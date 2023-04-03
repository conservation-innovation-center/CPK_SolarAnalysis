# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 16:06:03 2022

@author: mevans
"""

import geopandas as gpd
import requests
import pandas as pd
from os.path import join
import fiona 
from rasterstats import zonal_stats

network_drive = 'S:\\'
project_folder = 'SolarAnalysis\Data\Vector'

cpk_solar_file = join(network_drive, project_folder, 'CPK_arrays_GEEcovariates.geojson')
cpk_random_file = join(network_drive, project_folder, 'CPK_random_GEEcovariates.geojson')
cpk_states_file = join(network_drive, project_folder, 'CPK_states.shp')
roads_file= join(network_drive, project_folder, 'CPK_TIGER_roads.geojson')
parcel_file= join(network_drive, project_folder, 'parcels_CBW_bycounty.gdb')
tracts_file= join(network_drive, project_folder, 'CPK_TIGER_Tracts_2020_merge.shp')
padus_file= join(network_drive, project_folder, 'PADUS_Fee_CPK.shp')
lines_file = join(network_drive, project_folder, 'Electric_Power_Transmission_Lines/Electric_Power_Transmission_Lines.shp')
ssurgo_file = join(network_drive, 'SolarAnalysis\Data\Raster\SURRGO_4class.tif')

statesGDF = gpd.read_file(cpk_states_file).to_crs('epsg:3857')
arrayGDF = gpd.read_file(cpk_solar_file)
randomGDF = gpd.read_file(cpk_random_file)

# we can merge our arrays and random polygons
# randomGDF = randomGDF.rename(columns ={'label':'year'})
# arrayGDF = arrayGDF.rename(columns = {'year_left':'year'})
solarGDF = arrayGDF.append(randomGDF, ignore_index = True)
len(solarGDF) == len(randomGDF) + len(arrayGDF)
solarGDF = solarGDF.to_crs('epsg:3857')
# we will create an explicit 'index' column to track unique records
solarGDF['index'] = solarGDF.index

# we will use centroids for spatial join analyses to avoid problems with polygons overlapping multiple parcels
centroids = solarGDF.geometry.centroid
centroidGDF = gpd.GeoDataFrame(data = solarGDF['index'], geometry = centroids, crs = centroids.crs)

# we will use just polygons for near analyses
polysGDF = gpd.GeoDataFrame(data = solarGDF['index'], geometry = solarGDF.geometry, crs = solarGDF.crs)

# We'll deal with covariates one at a time
## census tracts
tractGDF = gpd.read_file(tracts_file)
# convert pop to mi-2
tractGDF['pdensity'] = tractGDF['POPULATION']/tractGDF['SQMI']

# join census tract FIPS and Population to arrays
# census tracts are from TIGER 2020 data
joined = centroidGDF.sjoin(tractGDF[['geometry', 'FIPS', 'POPULATION', 'pdensity']],
                        how = 'left',
                        predicate = 'within')

# all polygons should have a tract id - check for nas
# some will not, random rectangles landing in the water
# and polygons we detected in neighboring states due to chip buffer
joined.loc[joined['FIPS'].isnull(), ['index']]

# check for a 1 to 1 join
len(polysGDF) == len(joined)

# now join housing density data from census API
url = 'https://api.census.gov/data/2016/acs/acs5?get=B19001_001E,B00002_001E&for=tract:*&in=state:10,11,42,36,24,51'
  
income = 'B19001_001E'
housing_units = 'B00002_001E'
params = {'get': '{},{}'.format(income, housing_units), 'for': 'tract:*'}

response = requests.get(url = url)
response.url
dat = response.json()
df = pd.DataFrame(dat[1:], columns = dat[0])
df.columns = ['income', 'housing', 'state', 'county', 'tract']
df['FIPS'] = df.apply(lambda x: x['state']+x['county']+x['tract'], axis = 1)
df.head()

merged = joined.merge(df[['income', 'housing', 'FIPS']], on = 'FIPS', how = 'left')
len(joined) == len(merged)
# we don't need right index b/c we have FIPS as unique tract identifier
covGDF = merged.drop(['index_right'], axis =1)
# we only need to cary 'covGDF' forward
del(df, dat, tractGDF, joined, merged)

## PADUS
padusGDF = gpd.read_file(padus_file)

joined = covGDF.sjoin(
    padusGDF[['geometry', 'GAP_Sts']],
    how = 'left',
    predicate = 'within')

# check for 1 to 1
len(joined) == len(polysGDF)

# some arrays intersect multiple protected areas. we will take the most protective gap designation
dissolved = joined.dissolve(
    'index',
    aggfunc = {
        'GAP_Sts':'min',
        'FIPS':'first',
        'POPULATION':'first',
        'pdensity': 'first',
        'income':'first',
        'housing':'first',
        'index_right':'first'},
    as_index = False)

len(dissolved) == len(solarGDF)
covGDF = dissolved.rename({'index_right':'index_padus'}, axis = 1)
# we only need to cary dissolved forward
del(joined, dissolved, padusGDF)


## Parcels
# these came in a geodatabase with layers per county...gonna be weird

# get a list of layers in the gdb
counties = fiona.listlayers(parcel_file)


# loop through layers (i.e. counties) to get parcel size for each array polygon
GDFS= []
    
for county in counties[115:]:
    print(f'reading {county}')
    parcelGDF = gpd.read_file(parcel_file, driver = 'FileGDB', layer = county).to_crs(centroids.crs)
    if 'Shape_Area' not in parcelGDF.columns:
        parcelGDF['Shape_Area'] = parcelGDF.geometry.area
    joined = centroidGDF.sjoin(
        parcelGDF[['geometry', 'Shape_Area']],
        how = 'inner',
        predicate = 'within').reset_index()
    # add the county FIPS code for backwards traceability
    joined['county'] = county.split('_')[1]
    
    # we have to use this workaround to get the parcel id's corresponding to smallest parce
    dissolved = joined.loc[joined.groupby('index')['Shape_Area'].idxmin()]
    
    GDFS.append(dissolved)

parcelGDF = gpd.GeoDataFrame(pd.concat(GDFS, ignore_index = True)).rename({'index_right':'index_parcel'}, axis = 1)

# make sure nothing wierd happened while concatenating
len(parcelGDF) == sum([len(df) for df in GDFS])

# some parcels are in multiple county files, so we need to aggregate again
# here we will again use the smalles shape area 
dissolved = parcelGDF.loc[parcelGDF.groupby('index')['Shape_Area'].idxmin()]

merged = covGDF.merge(
    dissolved[['index', 'Shape_Area', 'county', 'index_parcel']],
    on = 'index',
    how = 'left')

if len(covGDF) == len(merged):
    covGDF = merged
else:
    print('merged geodataframe is of different length than polygon centroids')

del(merged, dissolved, parcelGDF, joined, GDFS)

## Latitude
covGDF['lat'] = centroids.y

## Transmission lines
linesGDF = gpd.read_file(lines_file)

nearest = polysGDF.sjoin_nearest(
    linesGDF[['geometry', 'VOLTAGE', 'ID']],
    how = 'left',
    max_distance = 75000,
    distance_col = 'line_dist').reset_index()

# check to make sure all polygons matched a power line. should be no nulls
len(nearest[nearest['index_right'].isnull()]) == 0

# some polygons match with two lines where they overlap. get the highest voltage
dissolved = nearest.loc[nearest.groupby('index')['VOLTAGE'].idxmax()]
# dissolved = nearest.dissolve('line_dist', aggfunc = 'min', as_index = False)

# ensure 1-to-1 join
if len(dissolved) == len(polysGDF):
    polysGDF = dissolved.drop(['level_0', 'index_right'], axis= 1)
else:
    print('joined geodataframe not same size as input')
    
# we only need to cary polysGDF forward
del(linesGDF, nearest, dissolved)

## Roads
# roads came from GEE with CRS:4326, need to reproject
roadGDF = gpd.read_file(roads_file).to_crs('epsg:3857')

nearest = polysGDF.sjoin_nearest(
    roadGDF[['geometry', 'mtfcc', 'linearid']],
    how = 'left',
    max_distance = 50000,
    distance_col = 'road_dist').reset_index()

# check to make sure all polygons matche a road. should be no nulls
len(nearest[nearest['index_right'].isnull()]) == 0

# check to see that we have a 1-to-1 join (we probably won't)
len(nearest) == len(polysGDF)

# some polygons match with two roads where they overlap pick the first
dissolved = nearest.loc[nearest.groupby('index')['index_right'].idxmin()]

if len(dissolved) == len(polysGDF):
    # we can drop index right b/c we have linearid from TIGER
    polysGDF = dissolved.drop(['index_right', 'level_0'], axis= 1)
else:
    print('joined geodataframe not same size as input')
    
# only need to carry merged forward
del(nearest, dissolved, roadGDF)

## Join the attributes we collected at centroids to those collected at polygons
vectorCovars = polysGDF.merge(
    covGDF.drop(['geometry'], axis = 1),
    on = 'index',
    how = 'inner')

# final checks
len(vectorCovars) == len(covGDF) == len(polysGDF)
len(vectorCovars[vectorCovars['index'].duplicated()]) == 0

# Finally, join our vector covariates back to our raster covariates
merged = solarGDF.merge(
    vectorCovars.drop(['geometry'], axis = 1),
    on = 'index',
    how = 'inner')

# final checks
len(merged) == len(solarGDF)
len(merged[merged['index'].duplicated()]) == 0


## SSURGO
# use zonal stats to extract the mean farmland importance value within each polygon

merged.to_file('./data/solar_analysis_data.geojson', driver = 'GeoJSON')

ssurgo_stats = zonal_stats(
    './data/solar_analysis_data.geojson',
    ssurgo_file,
    stats = 'mean')

#zonal_stats outputs a list of dictionaries, covert to array
ssurgo = [d['mean'] for d in ssurgo_stats]
merged['ssurgo'] = ssurgo

## Write our data out for analysis in R!
merged.drop(['geometry'], axis = 1).to_csv('./data/solar_analysis_data.csv')

no_income = merged[merged.income.isna()]
len(no_income)
fixed_income = []
for index, row in no_income.iterrows():
    neighbors = merged[merged.geometry.touches(row['geometry'])]
    income = neighbors.income[neighbors.income.notna()].astype(int).mean()
    FIPS = row.FIPS
    fixed_income.append(income)
    merged.loc[merged.FIPS == FIPS, 'income'] = income

## USGS GAP DATA
reptile_url = 'https://prod-is-s3-service.s3.amazonaws.com/ScienceBase/prod/5bf2eb54e4b045bfcae0c10b/a44a49612f4882091c4814bc48fbde3da0effa09/reptile_richness_habitat30m.tif?AWSAccessKeyId=AKIAI7K4IX6D4QLARINA&Expires=1675443607&Signature=x%2FvyFCPGVAPTvwE28XivOoQM0gw%3D'

reptile_stats = zonal_stats(
    solarGDF.to_crs(5070),
    reptile_url,
    band = 1,
    stats = 'mean')

reptile_list = [stat['mean'] for stat in reptile_stats]

taxaDF = pd.DataFrame(data = {'birds':bird_list, 'amphibians':amphibian_list, 'reptile':reptile_list, 'mammals':mammal_list})
taxaDF['state'] = solarGDF['state']
taxaDF['year'] = solarGDF['year']
taxaDF.to_csv('./data/biodiversity.csv')