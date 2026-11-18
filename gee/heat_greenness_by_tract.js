/**
 * Summer land surface temperature, greenness and land cover by census tract.
 *
 * Run in the Google Earth Engine Code Editor (https://code.earthengine.google.com).
 * Change the settings below, press Run, then start the export in the Tasks tab.
 * The CSV it writes is the input for scripts/run_analysis.py.
 *
 * Data
 *   Landsat 8 & 9 Collection 2 Level-2 (surface reflectance + surface temperature)
 *   ESA WorldCover 10 m v200 (tree cover, built-up, water)
 *   US Census TIGER 2020 tracts
 *
 * Author: Rajib Hossain
 */

// ---------------------------------------------------------------- settings
var STATE_FIPS = '13';            // Georgia
var COUNTY_FIPS = ['059'];        // Clarke County (Athens). Add more, e.g. ['121', '089'] for Atlanta
var START_YEAR = 2021;
var END_YEAR = 2025;
var SUMMER_MONTHS = [6, 8];       // June to August
var MAX_SCENE_CLOUD = 60;         // percent; pixels are cloud-masked anyway
var EXPORT_NAME = 'heat_greenness_tracts_GA_059';

// ---------------------------------------------------------------- study area
var tracts = ee.FeatureCollection('TIGER/2020/TRACT')
  .filter(ee.Filter.eq('STATEFP', STATE_FIPS))
  .filter(ee.Filter.inList('COUNTYFP', COUNTY_FIPS));
var region = tracts.geometry();

// ---------------------------------------------------------------- Landsat
function prepLandsat(img) {
  var qa = img.select('QA_PIXEL');
  var clear = qa.bitwiseAnd(1 << 1).eq(0)        // dilated cloud
    .and(qa.bitwiseAnd(1 << 3).eq(0))            // cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))            // cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0));           // snow
  var sr = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  var ndvi = sr.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
  var lst = img.select('ST_B10').multiply(0.00341802).add(149.0)
    .subtract(273.15).rename('LST_C');           // Kelvin -> Celsius
  return ndvi.addBands(lst).updateMask(clear)
    .copyProperties(img, ['system:time_start']);
}

var landsat = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .merge(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'))
  .filterBounds(region)
  .filter(ee.Filter.calendarRange(START_YEAR, END_YEAR, 'year'))
  .filter(ee.Filter.calendarRange(SUMMER_MONTHS[0], SUMMER_MONTHS[1], 'month'))
  .filter(ee.Filter.lt('CLOUD_COVER', MAX_SCENE_CLOUD))
  .map(prepLandsat);

print('Landsat scenes used', landsat.size());

// Median composite: robust to remaining cloud edges and one-off hot days
var summer = landsat.median();

// ---------------------------------------------------------------- land cover
var worldcover = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map');
var notWater = worldcover.neq(80);
var cover = ee.Image.cat([
  worldcover.eq(10).rename('tree_frac'),   // tree cover
  worldcover.eq(50).rename('built_frac')   // built-up
]);

// Exclude open water from temperature and NDVI (water is cool and would bias tracts with lakes)
var stack = summer.updateMask(notWater).addBands(cover);

// ---------------------------------------------------------------- zonal stats
var stats = stack.reduceRegions({
  collection: tracts,
  reducer: ee.Reducer.mean(),
  scale: 30,
  tileScale: 4
}).map(function (f) {
  var c = f.geometry().centroid(10).coordinates();
  return f.set({lon: c.get(0), lat: c.get(1)});
});

var columns = ['GEOID', 'NAME', 'ALAND', 'lon', 'lat', 'LST_C', 'NDVI', 'tree_frac', 'built_frac'];
Export.table.toDrive({
  collection: stats.select(columns, null, false),
  description: EXPORT_NAME,
  fileFormat: 'CSV',
  selectors: columns
});

// ---------------------------------------------------------------- map
Map.centerObject(region, 11);
var lstVis = {min: 28, max: 48, palette: ['313695', '74add1', 'ffffbf', 'f46d43', 'a50026']};
Map.addLayer(summer.select('LST_C').clip(region), lstVis, 'Summer LST (°C)');
Map.addLayer(summer.select('NDVI').clip(region), {min: 0, max: 0.8, palette: ['ffffff', '1a9850']}, 'NDVI', false);
Map.addLayer(ee.Image().paint(tracts, 0, 1), {palette: '000000'}, 'Tracts');

var tractLst = ee.Image().float().paint(stats, 'LST_C');
Map.addLayer(tractLst, lstVis, 'Tract mean LST (°C)', false);
