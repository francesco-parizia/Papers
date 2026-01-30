# MDSE extraction process ©SamueleDePetris ©FrancescoParizia

This code calculates iSCA and SP in order to obtain MDSE and nT
- 
This workflow is part of https://doi.org/10.5194/egusphere-2025-3028


-   [Earth Engine Homepage](https://earthengine.google.com/)
-   [Web Code Editor](https://code.earthengine.google.com/)
-   [GEE Catalog](https://developers.google.com/earth-engine/datasets/catalog/)



```javascript

var dem = ee.ImageCollection("COPERNICUS/DEM/GLO30"),
    modis = ee.ImageCollection("MODIS/061/MOD10A1"),
    geometry = /* color: #d63000 */ee.Geometry.Point([6.728941845258896, 45.120462038466634]);


/// Study area selection
var countries = ee.FeatureCollection("FAO/GAUL/2015/level1");
var Piemonte = countries.filter(ee.Filter.eq("ADM1_NAME", "Piemonte"));
var VdA = countries.filter(ee.Filter.eq("ADM1_NAME", "Valle D'aosta"));
var aoi = VdA.merge(Piemonte);

// 1. modis collection 

var mod10a1 = modis.filter(ee.Filter.date('2000-10-01', '2023-09-30')).filter(ee.Filter.bounds(aoi)); // NB: change years span <=========== (!)

// 2. downsampling dtm 
var dtm = dem.select('DEM').filter(ee.Filter.date('2000-10-01', '2023-09-30')).filter(ee.Filter.bounds(aoi)).mosaic().clip(aoi)
print(dtm,'dtm')
// Get the projection information for a band.
var b = mod10a1.select('NDSI_Snow_Cover').first();
var modisProjection = b.projection();
var dtm_res = dtm.reproject({crs: modisProjection, scale: 500});

// 3. Correcting area factor 

var slope = ee.Terrain.slope(dtm_res)
var caf = (ee.Image.constant(500*500).divide(slope.cos()).abs()).divide(10000) // area in ettari! 
print(caf)
Map.addLayer(caf,{min:0,max:1000},'caf',0);

// 4. Masking 


//create mask to extract only 'best' and 'good' quality data (values 0 and 1 in all the 16 bits)
//also create mask to extract pixels without inland water (bit0),
//no visible screen failure (bit1), no ndsi screen failure (bit2) and
//no solar zenith screen failure (bit7) 
//Refer:-   https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD10A1?hl=en
var filter = function(image){ 
  var basicQA = image.select('NDSI_Snow_Cover_Basic_QA');
  var basicBitMask = 1<<0|1<<1|1<<2|1<<3|1<<4|1<<5|1<<6|1<<7|1<<8|1<<9|1<<10|1<<11|1<<12|1<<13|1<<14|1<<15;
  var basicBitwiseResult = basicQA.bitwiseAnd(basicBitMask);
  var basicMask = basicBitwiseResult.eq(0).or(basicBitwiseResult.eq(1));
  var flagsQA = image.select('NDSI_Snow_Cover_Algorithm_Flags_QA');
  var inland = 1<<0;
  var visible = 1<<1;
  var ndsi = 1<<2;
  var solar = 1<<7;
  var Inland = (flagsQA.bitwiseAnd(inland)).eq(0);
  var Visible = (flagsQA.bitwiseAnd(visible)).eq(0);
  var NDSI = (flagsQA.bitwiseAnd(ndsi)).eq(0);
  var Solar = (flagsQA.bitwiseAnd(solar)).eq(0);
  var mask = basicMask.and(Inland).and(Visible).and(NDSI).and(Solar);
  image = image.updateMask(mask);
  return image;
};

//apply basic & flags QA mask to image collection
var mod_filtered = mod10a1.map(filter).select('NDSI_Snow_Cover');

print(mod_filtered,'MODIS collection flitered')



/// 5. gap-filling 

var mod_filtered = mod_filtered.map(function(image) {
  var timeImage = image.metadata('system:time_start').rename('timestamp')
  var timeImageMasked = timeImage.updateMask(image.mask().select(0))
  return image.addBands(timeImageMasked)
})



var days = 30
var millis = ee.Number(days).multiply(1000*60*60*24)

var maxDiffFilter = ee.Filter.maxDifference({
  difference: millis,
  leftField: 'system:time_start',
  rightField: 'system:time_start'
})

var lessEqFilter = ee.Filter.lessThanOrEquals({
  leftField: 'system:time_start',
  rightField: 'system:time_start'
})
var greaterEqFilter = ee.Filter.greaterThanOrEquals({
  leftField: 'system:time_start',
  rightField: 'system:time_start'
})

var filter1 = ee.Filter.and(maxDiffFilter, lessEqFilter)
 
var join1 = ee.Join.saveAll({
  matchesKey: 'after',
  ordering: 'system:time_start',
  ascending: false})
   
var join1Result = join1.apply({
  primary: mod_filtered,
  secondary: mod_filtered,
  condition: filter1
})

var filter2 = ee.Filter.and(maxDiffFilter, greaterEqFilter)
 
var join2 = ee.Join.saveAll({
  matchesKey: 'before',
  ordering: 'system:time_start',
  ascending: true})
   
var join2Result = join2.apply({
  primary: join1Result,
  secondary: join1Result,
  condition: filter2
})


var interpolateImages = function(image) {
  var image = ee.Image(image)
  // We get the list of before and after images from the image property
  // Mosaic the images so we a before and after image with the closest unmasked pixel
  var beforeImages = ee.List(image.get('before'))
  var beforeMosaic = ee.ImageCollection.fromImages(beforeImages).mosaic()
  var afterImages = ee.List(image.get('after'))
  var afterMosaic = ee.ImageCollection.fromImages(afterImages).mosaic()
  // Get image with before and after times
  var t1 = beforeMosaic.select('timestamp').rename('t1')
  var t2 = afterMosaic.select('timestamp').rename('t2')
  var t = image.metadata('system:time_start').rename('t')
  var timeImage = ee.Image.cat([t1, t2, t])
  var timeRatio = timeImage.expression('(t - t1) / (t2 - t1)', {
    't': timeImage.select('t'),
    't1': timeImage.select('t1'),
    't2': timeImage.select('t2'),
  })
  // Compute an image with the interpolated image y
  var interpolated = beforeMosaic
    .add((afterMosaic.subtract(beforeMosaic).multiply(timeRatio)))
  // Replace the masked pixels in the current image with the average value
  var result = image.unmask(interpolated)
  return result.copyProperties(image, ['system:time_start'])
}

var interpolatedCol = ee.ImageCollection(
  join2Result.map(interpolateImages))


// 6. Snow cover areas "real" stack 

var mult = function(image) {
  return image.multiply(caf).divide(100).set({'system:time_start': image.get('system:time_start')}).rename('sca');
};


var sca = (interpolatedCol.select('NDSI_Snow_Cover').map(mult)).select('sca')

// Define the chart and print it to the console. Achtung! no long series
var chart =
    ui.Chart.image
        .series({
          imageCollection: sca.select('sca'),
          region: geometry,
          reducer: ee.Reducer.mean(),
          scale: 500,
          xProperty: 'system:time_start'
        })
        .setSeriesNames(['sca'])
        .setOptions({
          title: 'SNOW cover area (ha) ',
          hAxis: {title: 'Date', titleTextStyle: {italic: false, bold: true}},
          vAxis: {
            title: 'SC area (ha)',
            titleTextStyle: {italic: false, bold: true}
          },
          lineWidth: 5,
          colors: ['e37d05', '1d6b99'],
          curveType: 'function'
        });
print(chart);

print(sca,'sca')
Map.addLayer(sca.median().clip(aoi),{min:0,max:50},'sca mediana',0);




// 7. integrating 



var dates = [970444800000,1001980800000] // set start-end in milliseconds of single year  <=========== (!)
var start = (ee.List(dates).get(0)) 
var end = (ee.List(dates).get(1))

var subset = sca.filterDate(start, end)

var sca_list = subset.sort('system:time_start').toList(subset.size())

var sca_integral = sca_list.slice(0,-1).zip(sca_list.slice(1))
  .map(function(f) { 
    var dy = ee.Image(ee.List(f).get(0))// get the left image  //ee.Image(ee.List(f).get(1)).subtract(
    var dx1 = ee.Number((ee.Image(ee.List(f).get(1))).get('system:time_start'))
    var dx2 = ee.Number((ee.Image(ee.List(f).get(0))).get('system:time_start'))
    var dx = ee.Image.constant(dx1.subtract(dx2)).divide(86400000)
    var d1 = (dy.multiply(dx)).abs()
    return(d1);
   })
var integral = ee.ImageCollection(sca_integral)
print(integral.first(),'integral first')

var sum = integral.sum().toInt32()



var palettes = require('users/gena/packages:palettes');
var palette = palettes.matplotlib.viridis[7];

Map.addLayer(sum.clip(aoi),{min:0,max:5000,palette: palette},'sum');




// 8. Snow duration 

var counting = function(image) {
  var background = ee.Image.constant(0)
  var mask = image.gt(background)
  var ii = image.updateMask(mask).multiply(0).add(1)
  var ss = ii.unmask()
  return ss.set({'system:time_start': image.get('system:time_start')}).rename('mask');
};

var binary = sca.map(counting).select('mask')

print(binary,'binary')



// doing using for loop for all years at same time  
var dates = [[970444800000,1001980800000],[1001980800000,1033516800000],[1033516800000,1065052800000],[1065052800000,1096588800000],[1096588800000,1128211200000],[1128211200000,1159747200000],[1159747200000,1191283200000],[1191283200000,1222819200000],[1222819200000,1254441600000],[1254441600000,1285977600000],[1285977600000,1317513600000],[1317513600000,1349049600000],[1349049600000,1380672000000],[1380672000000,1412208000000],[1412208000000,1443744000000],[1443744000000,1475280000000],[1475280000000,1506902400000],[1506902400000,1538438400000],[1538438400000,1569974400000],[1569974400000,1601510400000],[1601510400000,1633132800000],[1633132800000,1664668800000],[1664668800000,1696204800000]]


var fun = function(i) {
  var start = (ee.List(i).get(0))
  var end = (ee.List(i).get(1))
  var subset = binary.filterDate(start, end)

  var binary_list = subset.sort('system:time_start').toList(subset.size())

  var binary_integral = binary_list.slice(0,-1).zip(binary_list.slice(1))
    .map(function(f) { 
      var dy = ee.Image(ee.List(f).get(0))// get the left image  //ee.Image(ee.List(f).get(1)).subtract(
      var dx1 = ee.Number((ee.Image(ee.List(f).get(1))).get('system:time_start'))
      var dx2 = ee.Number((ee.Image(ee.List(f).get(0))).get('system:time_start'))
      var dx = ee.Image.constant(dx1.subtract(dx2).divide(86400000))
      var d1 = (dy.multiply(dx)).abs().set({'system:time_start': ee.Image(ee.List(f).get(0)).get('system:time_start')}).rename('mask1')
      return(d1);
     })
  var integral1 = ee.ImageCollection(binary_integral)
  var sum1 = integral1.sum()
  return sum1;
};

// Apply function to tile list
var filtered1 = dates.map(fun);
// Create image collection from image list
var filtered1 = ee.ImageCollection(filtered1);

print(filtered1,'fil')


Map.addLayer(filtered1.first().clip(aoi),{min:1000,max:10000},'integral first');


var stack_export1  = filtered1.toBands().toInt32();


Export.image.toDrive({
  image: stack_export1,
  folder: 'SnowCover',
  description: 'SP_Integral',     // <=========== (!)
  scale: 500, 
  maxPixels: 1e10,
  crs:'EPSG:32632',
  region: aoi
});


Export.image.toDrive({
  image: dtm_res,
  folder: 'SnowCover',
  description: 'DTM_Copernicus_500m',     // <=========== (!)
  scale: 500, 
  maxPixels: 1e10,
  crs:'EPSG:32632',
  region: aoi
});
