//Author : XU Jie, ITP CAS
//For details about this code, please refer to http://dx.doi.org/10.1080/01431161.2022.2118002
//Method: Jean-François Pekel1, Nature Letter 2016 and NIR+Otsu   
//Reference: https://gis.stackexchange.com/

var Hankou  = ee.FeatureCollection("users/christinaluintel/Hankou");

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//========================= 1.1 Preparison of SR data =========================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
var L8_SR = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
        .filterBounds(Hankou)
        .filterDate('2013-04-01','2019-12-22')
        .select(["SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","QA_PIXEL"],['blue','green','red','sr_nir','SWIR1','SWIR2','qa'])
        .map(function(image){
          var opticalBands = image.select('sr_nir','green','blue').multiply(0.0000275).add(-0.2);
          var date = image.date().format("YYYYMMdd");
          return image.addBands(opticalBands, null, true).set('date', date)
        });
        
//====================SR data quality ========================
//============================================================
var qualityfunction = function(image){
  var nir = image.select('sr_nir').gt(0)
  var blue = image.select('blue').gt(0)
  var red = image.select('red').gt(0)
  var green = image.select('green').gt(0)
  var SWIR1 = image.select('SWIR1').gt(0)
  var SWIR2 = image.select('SWIR2').gt(0)
  var quality = green.multiply(blue).multiply(red).multiply(nir).multiply(SWIR1).multiply(SWIR2).rename('quality')
  var qualitymin = quality.reduceRegion({
            reducer: ee.Reducer.min(),
            geometry: Hankou,
            scale: 30,
            });
  return image.set(qualitymin)
}
var L8_SR = L8_SR.map(qualityfunction).filter(ee.Filter.greaterThan('quality',0)).select('sr_nir','qa').sort('system:time_start');
print('L8_SR',L8_SR)

//=========== SR data cloud mask, avoid contamination of mean sir value ============
//=============== mean of nir under cloud would return null value ==================
//==================================================================================
var maskfunction = function(image){
  var image_qa = image.select('qa');
  // Create a mask for the dual QA bit "Cloud Confidence".
  // Bits 5-6: Cloud Confidence  
  // 0: Not Determined / Condition does not exist. 
  // 1: Low, (0-33 percent confidence)
  // 2: Medium, (34-66 percent confidence)
  // 3: High, (67-100 percent confidence)
  var RADIX = 2;  // Radix for binary (base 2) data.
  
  var extractQABits = function (qaBand, bitStart, bitEnd) {
    var numBits = bitEnd - bitStart + 1;
    var qaBits = qaBand.rightShift(bitStart).mod(Math.pow(RADIX, numBits));
    //Map.addLayer(qaBits, {min:0, max:(Math.pow(RADIX, numBits)-1)}, 'qaBits');
    return qaBits;
  };
  
  var bitStartCloudConfidence = 8;
  var bitEndCloudConfidence = 9;
  var qaBitsCloudConfidence = extractQABits(image_qa, bitStartCloudConfidence, bitEndCloudConfidence);
  // Test for clouds, based on the Cloud Confidence value.
  var mask_cloud_ori = qaBitsCloudConfidence.gte(2);
  var mask_cloud = ee.Image(0).select('constant').clip(Hankou).where(mask_cloud_ori.eq(1),1);
  // Get the internal_cloud_algorithm_flag bit.
  var sr_nir = image.select('sr_nir').updateMask(mask_cloud.not()).rename('sr_nir_cloudmasked')
  //mask those pixels from the image
  image = image.addBands(sr_nir);

  return image;
};
var L8_SR = L8_SR.map(maskfunction)

// https://developers.google.com/earth-engine/apidocs/ee-featurecollection-iterate
var L8_SR = L8_SR.select('sr_nir_cloudmasked').map(function(image) {
  var meanNir = image.select('sr_nir_cloudmasked').reduceRegion(
      {reducer: ee.Reducer.mean(), geometry: Hankou, scale: 30});
  return image.set(meanNir);
})
print(L8_SR)


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//========================= 1.2 Preparison of TOA ======================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
var L8_TOA = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')
        .filterBounds(Hankou)
        .filterDate('2013-04-01','2019-12-22')
        .select(["B2","B3","B4","B5","B6","B7","BQA"],['blue','green','red','nir','SWIR1','SWIR2','qa'])
        .map(function(image){
          var date = image.date().format("YYYYMMdd");
          return image.clip(Hankou).set('date', date);
        });

var qualityfunction = function(image){
  var nir = image.select('nir').gt(0)
  var blue = image.select('blue').gt(0)
  var red = image.select('red').gt(0)
  var green = image.select('green').gt(0)
  var quality = nir.multiply(blue).multiply(red).multiply(green)
  var qualitymin = quality.reduceRegion({
            reducer: ee.Reducer.min(),
            geometry: Hankou,
            scale: 30,
            });
  return image.set(qualitymin)
}
var L8_TOA = L8_TOA.map(qualityfunction).filter(ee.Filter.greaterThan('nir',0));
print('L8_TOA',L8_TOA)

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//=============== 1.3  Match the 'meanNir' Property with same date ==================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
var filter = ee.Filter.equals({
  leftField:  'date',
  rightField: 'date'
});
var simpleJoin = ee.Join.inner();
var innerJoin = ee.ImageCollection(simpleJoin.apply( L8_TOA,L8_SR, filter))
print('innerJoin',innerJoin)

function cleanJoin(feature){
  return ee.Feature(feature.get('primary')).copyProperties(feature.get('secondary'));
}
var L8 = innerJoin.map(cleanJoin);
print('L8',L8)

//--------------TOA ∩  SR =  AB∩CD, now L8 is BC, we also need A----------------
var simpleJoin = ee.Join.inverted();
var invertedJoin = ee.ImageCollection(simpleJoin.apply( L8_TOA,L8, filter))
print(invertedJoin)

L8 = L8.merge(invertedJoin)

var cloudfunction = function(image){
  //=====================================================================
  //--------------------mask out cloud ---------------------------------
  //get the Landsat 8 Pixel Quality Assessment(pixel_qa) Bit Index
  var cloud_qa = image.select('qa');
 
  // Create a mask for the dual QA bit "Cloud Confidence".
  // Bits 5-6: Cloud Confidence  
  // 0: Not Determined / Condition does not exist. \
  // 1: Low, (0-33 percent confidence)
  // 2: Medium, (34-66 percent confidence)
  // 3: High, (67-100 percent confidence)
  var RADIX = 2;  // Radix for binary (base 2) data.
  var extractQABits = function (qaBand, bitStart, bitEnd) {
    var numBits = bitEnd - bitStart + 1;
    var qaBits = qaBand.rightShift(bitStart).mod(Math.pow(RADIX, numBits));
    //Map.addLayer(qaBits, {min:0, max:(Math.pow(RADIX, numBits)-1)}, 'qaBits');
    return qaBits;
  };
  
  var bitStartCloudConfidence = 5;
  var bitEndCloudConfidence = 6;
  var qaBitsCloudConfidence = extractQABits(cloud_qa, bitStartCloudConfidence, bitEndCloudConfidence);
  // Test for clouds, based on the Cloud Confidence value.
  var mask_cloud_ori = qaBitsCloudConfidence.gte(2)
  var mask_cloud     = ee.Image(0).select('constant').clip(Hankou).where(mask_cloud_ori.eq(1),1).rename('mask_cloud').reproject({crs:'EPSG:2163',scale:30});
  //mask those pixels from the image
  
  bitStartCloudConfidence = 14;
  bitEndCloudConfidence = 15;
  qaBitsCloudConfidence = extractQABits(cloud_qa, bitStartCloudConfidence, bitEndCloudConfidence);
  // Test for cirrus, based on the Cirrus Confidence value.
  var mask_cirrus_ori = qaBitsCloudConfidence.gte(2)
  var mask_cirrus     = ee.Image(0).select('constant').clip(Hankou).where(mask_cirrus_ori.eq(1),1).rename('mask_cirrus').reproject({crs:'EPSG:2163',scale:30});
  
  //-------------------- get the cirrus area ---------------------------------
  var area = ee.Image.pixelArea().reproject({crs:'EPSG:2163',scale:30}); 
  mask_cirrus = mask_cirrus.eq(1);
  var cirrusArea = mask_cirrus.multiply(area).rename('cirrusArea');
  
  var stats = cirrusArea.reduceRegion({
    //crs: 'EPSG:4326',
    reducer: ee.Reducer.sum(), 
    geometry: Hankou, 
    scale: 30,
  });
  
  image = image.addBands(mask_cloud).set('cirrusratio', ee.Number(stats.get('cirrusArea')).divide(ee.Number(Hankou.geometry().area())))
  return image;
};
var L8 = L8.map(cloudfunction);
print('L8',L8);


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//============================= 2.1 get Pekel data ==================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
var L8_Pekel = L8.filter(ee.Filter.calendarRange(5,9,'month'))
                .merge(L8.filter(ee.Filter.calendarRange(10,12,'month')).merge(L8.filter(ee.Filter.calendarRange(1,4,'month'))).filter(ee.Filter.greaterThan('sr_nir_cloudmasked',0.115)))
                .merge(L8.filter(ee.Filter.calendarRange(10,12,'month')).merge(L8.filter(ee.Filter.calendarRange(1,4,'month'))).filter(ee.Filter.lessThan('sr_nir_cloudmasked',0.115)).filter(ee.Filter.greaterThan('cirrusratio',0.2)))
                .sort('system:time_start')  // Sort chronologically in ascending order.
print('L8_Pekel',L8_Pekel)
// parameter function:
var maskfunction_Pekel = function(image){
  //=================================================================
  //------------1. mask out water/cloud pixel -----------------------
  var NDVI = image.normalizedDifference(['nir','red']).rename('nd');
  var MNDWI = image.normalizedDifference(['green','nir']).rename('NDWI');
  var hsv = image.select(['SWIR2', 'nir', 'red']).rgbToHsv();   // Convert the RGB bands to the HSV color space.
  image = image.addBands(NDVI).addBands(MNDWI).addBands(hsv)
  var mask_water_ori = image.expression(
      "(b('nd') < 0&&b('value') < 0.62&&(((b('hue')<((-9.867784585617413*b('nd'))+238.26034242940045))"+
      "&&(b('hue')>((-12960.000000000335*b('nd'))-12714.048607819708))||(b('hue')>((23.627546071775214*b('nd'))+255.53176874753507)))"+
      "||((b('hue')<((-54.685799109352004*b('nd'))+215.15052322834936))&&(b('hue')<((23.627546071775214*b('nd'))+255.53176874753507))"+
      "&&(b('hue')>((-7.321079389910027*b('nd'))+224.6166270396205)))||((b('hue')<((-172.0408163265306*b('nd'))+191.69646750224035))"+
      "&&(b('hue')<((-7.321079389910027*b('nd'))+224.6166270396205))&&(b('hue')>((-38.11764705882351*b('nd'))+193.8533786110101)))"+
      "||((b('hue')>((-52.06378986866776*b('nd'))+179.92232432949075))&&(b('hue')<((-879.6226415094455*b('nd'))+180.3004476242325))"+
      "&&(b('hue')<((-38.11764705882351*b('nd'))+193.8533786110101)))))? 1" + 
      //"(b('hue')>((23.627546071775214*b('nd'))+255.53176874753507))?1"+ this condition makes the image full balck.
      // So change && to || in front of it.

      ":(b('nd') > 0&&b('value')<0.62&&(((b('hue')<((-119.15098406819945*b('nd'))+180.0533162435398))"+
      "&&(b('hue')>((-994.2857142867327*b('nd'))+180.04805813312743))&&(b('hue')>((-116.5000234173271*b('nd'))+179.9633248496054)))"+
      "||((b('hue')<((-2368.4258422651174*b('nd'))+256.40879883589054))&&(b('hue')<((-116.5000234173271*b('nd'))+179.9633248496054))"+
      "&&(b('hue')>((-267.6720052547653*b('nd'))+179.97791758964533)))"+
      "||((b('hue')<((-108.07947019867622*b('nd'))+179.67747476669464))&&(b('hue')>((-2368.4258422651174*b('nd'))+256.40879883589054))"+
      "&&(b('hue')>((58.99660016815455*b('nd'))+168.09286521078695)))||((b('hue')<((-104.45621862799788*b('nd'))+179.4262481567021))"+
      "&&(b('hue')<((58.99660016815455*b('nd'))+168.09286521078695))&&(b('hue')>((-52.1565190088889*b('nd'))+172.13690440390852)))"+
      "||((b('hue')<((-52.1565190088889*b('nd'))+172.13690440390852))&&(b('hue')>((-204.2258047185466*b('nd'))+177.66958001421082))"+
      "&&(b('hue')>((37.74894387447151*b('nd'))+159.60620482085795)))))? 1" + //(b('nd') < 0&&

      ":(b('SWIR1')<0||b('SWIR2')<0)? 1"+
      ": 0")
  var mask_water = ee.Image(0).clip(Hankou).where(mask_water_ori.eq(1),1).rename('mask_water').reproject({crs:'EPSG:2163',scale:30}) 

  //-------------------- get the 900m2 area band---------------------------------
  var area = ee.Image.pixelArea().reproject({crs:'EPSG:2163',scale:30}); 
  var water01 = mask_water.eq(1);
  var oriArea = water01.multiply(area).rename('oriArea');
  
  var stats = oriArea.reduceRegion({
    //crs: 'EPSG:4326',
    reducer: ee.Reducer.sum(), 
    geometry: Hankou, 
    scale: 30,
  });

  //mask those pixels from the image
  image = image.addBands(mask_water).addBands(oriArea).set(stats);

  return image;
};

var L8_mask_Pekel = L8_Pekel.map(maskfunction_Pekel);
print('L8_mask_Pekel',L8_mask_Pekel);


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//============================= 2.1 get NIR data ==================================
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
var L8_NIR = ee.ImageCollection(simpleJoin.apply( L8,L8_Pekel, filter)).sort('system:time_start')  // Sort chronologically in ascending order.
print('L8_NIR',L8_NIR)
var histofunction = function(image){
  // Compute the histogram of the NIR band.  The mean and variance are only FYI.
  var histogram = image.select('nir').reduceRegion({
    reducer: ee.Reducer.histogram(255, 2)
      .combine('mean', null, true)
      .combine('variance', null, true), 
    geometry: Hankou, 
    scale: 30,
    bestEffort: true
  });
  image = image.set(histogram)
  return image;
}
var L8_histo = L8_NIR.map(histofunction);
print('L8_histo',L8_histo)

var maskfunction_NIR = function(image){
  // Return the DN that maximizes interclass variance in B5 (in the region).
  var otsu = function(histogram) {
    var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
    var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
    var size = means.length().get([0]);
    var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
    var mean = sum.divide(total);
    var indices = ee.List.sequence(1, size);
  
    // Compute between sum of squares, where each mean partitions the data.
    var bss = indices.map(function(i) {
      var aCounts = counts.slice(0, 0, i);
      var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
      var aMeans = means.slice(0, 0, i);
      var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
      var bCount = total.subtract(aCount);
      var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
      return aCount.multiply(aMean.subtract(mean).pow(2)).add(
          bCount.multiply(bMean.subtract(mean).pow(2)));
    });
    // Return the mean value corresponding to the maximum BSS.
    return means.sort(bss).get([-1]);
  };

  var threshold = otsu(image.get('nir_histogram'));
  var mask_water_ori = image.select('nir').lt(threshold)
  var mask_water = ee.Image(0).clip(Hankou).where(mask_water_ori.eq(1),1).rename('mask_water').reproject({crs:'EPSG:2163',scale:30}) 

  //-------------------- get the 900m2 area band---------------------------------
  var area = ee.Image.pixelArea().reproject({crs:'EPSG:2163',scale:30}); 
  var water01 = mask_water.eq(1);
  var oriArea = water01.multiply(area).rename('oriArea');
  
  var stats = oriArea.reduceRegion({
    //crs: 'EPSG:4326',
    reducer: ee.Reducer.sum(), 
    geometry: Hankou, 
    scale: 30,
  });

  //mask those pixels from the image
  image = image.addBands(mask_water).addBands(oriArea).set(stats);

  return image;
};

var L8_mask_NIR = L8_histo.filter(ee.Filter.gt('nir_mean', 0)).map(maskfunction_NIR);
print('L8_mask_NIR',L8_mask_NIR);

// merge collections
var L8_mask = L8_mask_Pekel.merge(L8_mask_NIR).sort("system:time_start")
print('L8_mask',L8_mask)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//------------  2. reduce the image collection to one image by summing all the rasters  ---------
//--------------------   and select out the minum count in water area    ------------------
var Hankou_water  = ee.FeatureCollection("users/christinaluintel/Hankou100751water");
var Pekel_pixelcounts = ee.Image(L8_mask.filter(ee.Filter.calendarRange(5,9,'month')).select('mask_water').sum().rename('counts'))
                        .reproject({crs:'EPSG:2163',scale:30}).toInt()
var Pekel_min =  Pekel_pixelcounts                       
                        .reduceRegion({
                            reducer: ee.Reducer.min(),
                            geometry: Hankou_water,
                            scale: 30,
                        });
var NIR_pixelcounts = ee.Image(L8_mask.filter(ee.Filter.calendarRange(10,12,'month')).merge(L8_mask.filter(ee.Filter.calendarRange(1,4,'month'))).select('mask_water').sum().rename('counts'))
                        .reproject({crs:'EPSG:2163',scale:30}).toInt()
var NIR_min =  NIR_pixelcounts                        
                        .reduceRegion({
                            reducer: ee.Reducer.min(),
                            geometry: Hankou_water,
                            scale: 30,
                        });

var L8_mask_Pekel = L8_mask.filter(ee.Filter.calendarRange(5,9,'month')).map(function(image){
  var Pekel_threshold = ee.Image(ee.Number(Pekel_min.get('counts'))).select('constant').rename('threshold').clip(Hankou).reproject({crs:'EPSG:2163',scale:30});
  return image.addBands(Pekel_pixelcounts).addBands(Pekel_threshold)
})

var L8_mask_NIR = L8_mask.filter(ee.Filter.calendarRange(10,12,'month')).merge(L8_mask.filter(ee.Filter.calendarRange(1,4,'month'))).map(function(image){
  var NIR_threshold = ee.Image(ee.Number(NIR_min.get('counts'))).select('constant').rename('threshold').clip(Hankou).reproject({crs:'EPSG:2163',scale:30});
  return image.addBands(NIR_pixelcounts).addBands(NIR_threshold)
})

var L8_mask = L8_mask_Pekel.merge(L8_mask_NIR).sort("system:time_start")
print('L8_mask',L8_mask)

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//--------------- 3. fill the water pixels under cloud and boats-----------------------------------
var fillfunction = function(image){
  var counts = image.select('counts')
  var threshold = image.select('threshold')
  var green = image.select('nir')
  var stripe = green.eq(-9999).updateMask(green).rename('stripe')   //null value in stripe part, main parts becomes zero
  var backgroud1 = ee.Image(0).select('constant').clip(Hankou)
  var vispart    = backgroud1.where(stripe.eq(0),1)
  
  var mask_water = image.select('mask_water')
  var mask_cloud = image.select('mask_cloud')
  var counts3 = ee.Image(0).clip(Hankou).rename('counts3')
  var counts1 = counts3.where(stripe.eq(0),counts).rename('counts1') // main parts have sum values, null can't select out.
  var counts2 = counts3.where(counts1.eq(0),counts).rename('counts2')  // stripe parts have sum values
  var mask_water_bg  = counts3.where(mask_water.eq(1),1).rename('mask_water_bg')
  counts3 = counts3.addBands(counts1).addBands(mask_cloud).addBands(mask_water).addBands(counts2).addBands(mask_water_bg).addBands(threshold)
  var mask_water_main    = counts3.expression("(b('counts1')>b('threshold')&&b('mask_cloud')>0&&b('mask_water')<1)?1"+
                                          ": b('mask_water_bg')").rename('mask_water_main').reproject({crs:'EPSG:2163',scale:30})   //main parts gap filling  
  counts3 = counts3.addBands(mask_water_main)
  //var mask_water_whole   = mask_water_main.where((mask_water_main.eq(0)&&counts2.gt(120)),1).rename('mask_water_whole') //stripe parts gap filling  // it seems that where have some bug that doesn't work well, I don't know why.
  var mask_water_whole = counts3.expression("(b('mask_water_main')<1&&b('counts2')>b('threshold'))?1"+ 
                                            ": b('mask_water_main')").rename('mask_water_whole')  //stripe parts gap filling
  
  // residual non-water pixels 
  var mask = mask_water_whole.eq(0)
  var masked = mask.updateMask(mask).rename('masked')
  var patchsize = masked.connectedPixelCount({
    maxSize: 25, eightConnected: false
  }).rename('patchsize');  
  mask_water_whole = mask_water_whole.addBands(patchsize)
  var fill_mask = mask_water_whole.expression("((b('patchsize')<25&&b('mask_water_whole')<1)?1"+ 
                                            ": b('mask_water_whole'))").reproject({crs:'EPSG:2163',scale:30}).rename('fill_mask')  //block parts gap filling

  //-------------------- get visArea ---------------------------------
  var area = ee.Image.pixelArea().reproject({crs:'EPSG:2163',scale:30}); 
  vispart = vispart.eq(1);
  var visArea = vispart.multiply(area).rename('visArea');
  var vis = visArea.reduceRegion({
    //crs: 'EPSG:4326',
    reducer: ee.Reducer.sum(), 
    geometry: Hankou, 
    scale: 30,
  });
  
  //-------------------- get water area in vis part --------------------------------
  var water02 = mask_water_main.eq(1);
  var mainArea = water02.multiply(area).rename('mainArea');  
  var main = mainArea.reduceRegion({
    //crs: 'EPSG:4326',
    reducer: ee.Reducer.sum(), 
    geometry: Hankou, 
    scale: 30,
  });
  
  //-------------------- get filled Area -------------------------
  var water03 = fill_mask.eq(1);
  var fillArea = water03.multiply(area).rename('fillArea');
  var fill = fillArea.reduceRegion({
    reducer: ee.Reducer.sum(), 
    geometry: Hankou, 
    scale: 30,
  });
  
  //---------------- get permanent Area -------------------
  var permanent_water = counts.gt(image.select('threshold'));
  var permanentArea = permanent_water.multiply(area).rename('permanentArea')
  var permanent = permanentArea.reduceRegion({
    //crs: 'EPSG:4326',
    reducer: ee.Reducer.sum(), 
    geometry: Hankou, 
    scale: 30,
  });
  
  return image.addBands(fill_mask).addBands(fillArea).set(main).set(vis).set(permanent).set(fill);
}
var L8_fill = L8_mask.map(fillfunction)
print('L8_fill',L8_fill)

//======================6. overlayed/permanent>0.9, oriArea/fillArea>0.9 =======================
//-------------------  caculate the filled and permanent overlayed area --------------------------
var overlayfunction = function(image){
  var overlay = image.expression("(b('counts')>b('threshold')&&b('fill_mask')>0?1"+
                                  ": 0)").rename('overlay').reproject({crs:'EPSG:2163',scale:30});
  var area = ee.Image.pixelArea().reproject({crs:'EPSG:2163',scale:30});                                
  var overlay_pixel = overlay.eq(1);
  var overlayArea = overlay_pixel.multiply(area).rename('overlayArea');
  var overlaystats = overlayArea.reduceRegion({
    reducer: ee.Reducer.sum(), 
    geometry: Hankou, 
    scale: 30,
  });
  return image.set(overlaystats); 
}
var L8_overlay = L8_fill.map(overlayfunction) 

var ratiofunction = function(image){
  return image.setMulti({
    overlayratio: ee.Number(image.get('overlayArea')).divide(ee.Number(image.get('permanentArea'))), 
    orifullratio: ee.Number(image.get('oriArea')).divide(ee.Number(image.get('mainArea'))), 
    ratio: ee.Number(image.get('visArea')).divide(ee.Number(Hankou.geometry().area()))
  })
}
var L8_overlay = L8_overlay.map(ratiofunction)
print('L8_overlay',L8_overlay)

//--------------- 2. filter out those with no water pixel detection -------------------------
var L8_water  = L8_overlay.filter(ee.Filter.greaterThan('fillArea',900)) //water area with greater than one pixel
                          .filter(ee.Filter.greaterThan('orifullratio',0.9)) //original-water/fill-water >0.9
                          .filter(ee.Filter.greaterThan('overlayratio',0.9)) //fill-water/permanent-water >0.9 
                          .select('fill_mask','fillArea')

print('L8_water',L8_water)

// //==========================================================================================

//--------------- 3.2 add the maxium area to the image collection------------
//                    Convert the zones of water to vectors.
  var maxareafunction = function(image){
    var vector = image.reduceToVectors({
      geometry: Hankou,
      scale: 1,
      geometryType: 'polygon',
      eightConnected: false,
      labelProperty: 'zone',
      reducer: ee.Reducer.sum()
    });
    var maxwater = vector.aggregate_array('sum')
    var listb = ee.List(maxwater)
    var list  = listb.sort()
    var counts = vector.size().subtract(1)
    var highest =  list.get(counts);
    return image.set('maxwaterArea', highest);
  }
    
  var L8_maxarea = L8_water.map(maxareafunction);
  print('L8_maxarea',L8_maxarea)
  
  var ratio2function = function(image){
  return image.setMulti({
    maxratio: ee.Number(image.get('maxwaterArea')).divide(ee.Number(image.get('permanentArea'))), 
  })
}
var L8_final = L8_maxarea.map(ratio2function)
var L8_final  = L8_final.filter(ee.Filter.greaterThan('maxratio',900)) //maxwater/permanent-water >0.9 
print('L8_final',L8_final)

Export.table.toDrive({
  collection: L8_maxarea,
  description: 'L8_100715_ALL',
  folder: 'GEE_geohackweek',
  fileFormat: 'CSV'
});
