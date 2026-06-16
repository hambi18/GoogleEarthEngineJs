Map.addLayer(aoi);
Map.centerObject(aoi, 9);

/////////////////////////YEAR SEQUENCE////Ikongo sorona///////////////////////////////////////
// Define start and end years.
var startYear = 1990;
var endYear = 2020;

/////////////////////////////////////////////////////////////////////////////////////////////////
// --------------------------------CLOUD MASK FOR COLLECTION 2 TIER 1-------------------------------------
// QA_PIXEL bit flags for Collection 2
// Bit 3: Cloud
// Bit 4: Cloud Shadow
// Bit 1: Dilated Cloud
// Bit 2: Cirrus Cloud

function maskCloudsC2(img) {
  var qa = img.select('QA_PIXEL');
  var dilatCloudBitMask = (1 << 1);      // Dilated Cloud
  var cirrusCloudBitMask = (1 << 2);     // Cirrus Cloud
  var cloudBitMask = (1 << 3);           // Cloud
  var cloudShadowBitMask = (1 << 4);     // Cloud Shadow
  
  var mask = qa.bitwiseAnd(dilatCloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusCloudBitMask).eq(0))
    .and(qa.bitwiseAnd(cloudBitMask).eq(0))
    .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0));
  
  return img.updateMask(mask);
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Define coefficients supplied by Roy et al. (2016) for translating ETM+
// surface reflectance to OLI surface reflectance.
var coefficients = {
  itcps: ee.Image.constant([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172]).multiply(10000),
  slopes: ee.Image.constant([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071])
};

//-----------------------------------------STANDARDIZATION - OLI- ETM+- TM---------------------------------------

// 1_Define function to get and rename bands of interest from OLI/L8 Collection 2 Tier 1.
function renameOLI(img) {
  return img.select(
    ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
  ).multiply(0.0000275).add(-0.2);
}

// 2_Define function to get and rename bands of interest from ETM+/L7 Collection 2 Tier 1.
function renameETM(img) {
  return img.select(
    ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
  ).multiply(0.0000275).add(-0.2);
}

// 3_Define function to get and rename bands of interest from TM/L5 Collection 2 Tier 1.
function renameTM(img) {
  return img.select(
    ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
  ).multiply(0.0000275).add(-0.2);
}

//------------------------------- HARMONIZATION TRANSFORMATION------------------------------------------
// Define function to apply harmonization transformation.
// L7 to L8
function etm2oli(img) {
  return img.select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
    .multiply(coefficients.slopes)
    .add(coefficients.itcps)
    .round()
    .toShort();
}

// L5 to L8
function tm2oli(img) {
  return img.select(['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'])
    .multiply(coefficients.slopes)
    .add(coefficients.itcps)
    .round()
    .toShort();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//------------------------IMAGE PREPARATION FUNCTION OLI-ETM+-TM---------------------------------
// Define function to prepare OLI or L8 Collection 2 Tier 1 images.
function prepOLI(img) {
  var orig = img;
  img = maskCloudsC2(img);
  img = renameOLI(img);
  return ee.Image(img).copyProperties(orig, orig.propertyNames());
}

// Define function to prepare ETM+ or L7 Collection 2 Tier 1 images.
function prepETM(img) {
  var orig = img;
  img = maskCloudsC2(img);
  img = renameETM(img);
  img = etm2oli(img);
  return ee.Image(img).copyProperties(orig, orig.propertyNames());
}

// Define function to prepare TM or L5 Collection 2 Tier 1 images.
function prepTM(img) {
  var orig = img;
  img = maskCloudsC2(img);
  img = renameTM(img);
  img = tm2oli(img);
  return ee.Image(img).copyProperties(orig, orig.propertyNames());
}

// Gap fill function for Landsat 7 (SLC-off)
function gapFill(image) {
  var filled = image.focal_mean(2, 'square', 'pixels', 1);
  return filled.blend(image);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////IMAGE COLLECTION//////////////////////////////////////////////////////////
// Updated to use Collection 2 Tier 1 - T1 indicates Tier 1 data

// Landsat 5 Surface Reflectance Collection 2 Tier 1
var tmCol = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterDate('1990-01-01', '2011-12-31')
  .filter(ee.Filter.lt('CLOUD_COVER', 10))
  .filterBounds(aoi)
  .map(function(image) {
    return image.clip(aoi);
  })
  .map(function(a) {
    return a.set('year', ee.Image(a).date().get('year'));
  });

// Landsat 7 Surface Reflectance Collection 2 Tier 1
var etmCol = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterDate('2012-01-01', '2012-12-31')
  .filter(ee.Filter.lt('CLOUD_COVER', 10))
  .filterBounds(aoi)
  .map(function(image) {
    return image.clip(aoi);
  })
  .map(function(a) {
    return a.set('year', ee.Image(a).date().get('year'));
  })
  .map(gapFill); // Gap fill for SLC-off

// Landsat 8 Surface Reflectance Collection 2 Tier 1
var oliCol = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterDate('2013-01-01', '2020-12-31')
  .filter(ee.Filter.lt('CLOUD_COVER', 10))
  .filterBounds(aoi)
  .map(function(image) {
    return image.clip(aoi);
  })
  .map(function(a) {
    return a.set('year', ee.Image(a).date().get('year'));
  });

print('tmCol', tmCol.size());
print('etmCol', etmCol.size());
print('oliCol', oliCol.size());

// Collection filter
var colFilter = ee.Filter.and(
  ee.Filter.bounds(aoi),
  ee.Filter.calendarRange(1, 365, 'day_of_year'),
  ee.Filter.lt('GEOMETRIC_RMSE_MODEL', 10)
);

//---------------------------MERGING OLI-ETM+-TM-------------------------------------
// Filter collections and prepare them for merging
oliCol = oliCol.filter(colFilter).map(prepOLI);
etmCol = etmCol.filter(colFilter).map(prepETM);
tmCol = tmCol.filter(colFilter).map(prepTM);

// MERGE L5, L7 & L8
var collection_merge = ee.ImageCollection(tmCol.merge(etmCol.merge(oliCol)));
print('collection_full', collection_merge);
print('collection_full size', collection_merge.size());

///USER INPUT///////
// Create a list of years to be iterated over
var years = ee.List.sequence(1990, 2020);
print('list of years', years);

// Create a collection with 1 image for each year (median composite)
var collectYear = ee.ImageCollection(years
  .map(function(y) {
    var start = ee.Date.fromYMD(y, 1, 1);
    var end = start.advance(12, 'month');
    return collection_merge.filterDate(start, end).reduce(ee.Reducer.median());
  }));
print('col with 1 image for each year', collectYear);

// Count number of bands in each image, if 0 remove from image collection
var nullimages = collectYear
  .map(function(image) {
    return image.set('count', image.bandNames().length());
  })
  .filter(ee.Filter.gt('count', 0));
print('col with 1 img for each year with bands', nullimages);

// Get image for display and training
var median_image = ee.Image(nullimages.first());
print('median_image info', median_image);
print('median_image bands', median_image.bandNames());

Map.addLayer(median_image, {
  bands: ['Red_median', 'Green_median', 'Blue_median'],
  min: 0,
  max: 0.3
}, 'Median Composite');

Map.centerObject(aoi, 9);
Map.addLayer(aoi, {color: 'yellow'}, 'AOI');

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //------------------------------CLASSIFICATION --------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Merge training points/polygons (assuming Forest and Nonforest are pre-defined geometries)
var training_points = Forest.merge(Nonforest);
print('Training points count:', training_points.size());

// Extract training data from median image
var training_data = median_image.sampleRegions({
  collection: training_points,
  properties: ['LC'],
  scale: 30
});

print('Training data', training_data);
print('Training data size', training_data.size());
print('First training point:', training_data.first());

// Build a random forest classifier with 10 trees
var trainingclassifier = ee.Classifier.smileRandomForest(10)
  .train({
    features: training_data,
    classProperty: 'LC',
    inputProperties: ['Red_median', 'Green_median', 'Blue_median', 'NIR_median', 'SWIR1_median', 'SWIR2_median']
  });

print('Training classifier:', trainingclassifier);

// Print Confusion Matrix and Overall Accuracy
var confusionMatrix = trainingclassifier.confusionMatrix();
print('Error matrix: ', confusionMatrix);
print('Training overall accuracy: ', confusionMatrix.accuracy());

// Apply the trained classifier to the median image
var classified = median_image.select(['Red_median', 'Green_median', 'Blue_median', 'NIR_median', 'SWIR1_median', 'SWIR2_median'])
  .classify(trainingclassifier);

// Visualization parameters for classification
var classVisParam = {
  min: 0,
  max: 1,
  palette: ['006400', 'ff0000'] // Dark green for Forest, Red for Non-forest
};

Map.addLayer(classified, classVisParam, "Classified Map");

print('Classification complete!');
