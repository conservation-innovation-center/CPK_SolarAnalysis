# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 11:47:09 2022

@author: mevans
"""

import geopandas as gpd
import pandas as pd
from os.path import join
from os import listdir, walk
import fiona 

network_drive = 'K:\\'
project_folder = 'GIS\CBP_Obj_1\data\\20172018_County_Data\County_Planimetrics_Projected'

# get the list of counties in the gdb that we previously processed
parcel_file= join(network_drive, project_folder, 'parcels_CBW_bycounty.gdb')
counties = fiona.listlayers(parcel_file)


merged = gpd.read_file('./data/solar_analysis_data.geojson', driver = 'GeoJSON')
# just focus on records for which we don't have parcel data yet
missing = merged[merged['Shape_Area'].isnull()]
centroids = missing.geometry.centroid
centroidGDF = gpd.GeoDataFrame(data = missing['index'], geometry = centroids, crs = centroids.crs)

def get_parcel_area(parcelGDF, file):
    if 'Shape_Area' not in parcelGDF.columns:
        parcelGDF['Shape_Area'] = parcelGDF.geometry.area
    joined = centroidGDF.sjoin(
        parcelGDF[['geometry', 'Shape_Area']],
        how = 'inner',
        predicate = 'within').reset_index()
    # add the county FIPS code for backwards traceability
    joined['county'] = file.split('_')[1]
    
    # we have to use this workaround to get the parcel id's corresponding to smallest parce
    dissolved = joined.loc[joined.groupby('index')['Shape_Area'].idxmin()]
    return(dissolved)
    
# create a dictionary to store new gdf per county
dic = {}

# MARYLAND


# first we'll get counties not in the gdb
# rstrip removes trailing '_parcels' to our names match files in md directory
countylist = [county.rstrip('_parcels') for county in counties if '24' in county]
countylist.sort()

# list of county folder names
files = listdir(join(network_drive, project_folder, 'MD'))
# return files that are not in the list of those we've dealt with already
filesToProcess = [file for file in files + countylist if file not in countylist]
filesToProcess

shapefiles = []
for root, dirs, files in walk(join(network_drive, project_folder, 'MD')):
    # if we have reached files
    if len(files) > 0:
        # find files with 'parcel' and '.shp'
        file = [file for file in files if 'arcel' in file and '.shp' in file]
        if len(file) > 0:
            shapefiles.append(join(root, file[0]))

finalshapefiles = {}
for county in filesToProcess:
    for file in shapefiles:
        if county in file:
            finalshapefiles[county]=file

processedFiles = []

for county, file in finalshapefiles.items():
    print(county)
    print(file)
    parcelGDF = gpd.read_file(file).to_crs(merged.crs)
    temp = get_parcel_area(parcelGDF, county)
    dic[county] = temp
    processedFiles.append(county)
    
# we wont find shapefiles that are part of gdb, so do these individually    
[file for file in filesToProcess if file not in processedFiles]

# NOTE: NO PARCEL DATA FOR ceci_24015??

## DELAWARE
# first we'll get counties not in the gdb
# rstrip removes trailing '_parcels' to our names match files in md directory
decountylist = [county.rstrip('_parcels') for county in counties if '10' in county]
decountylist.sort()

# we only missed one county in DE
file = join(network_drive, project_folder, 'DE', 'suss_10005\CIC_Data\CIC_Data\parcels_proj_MPtoSP.shp')

parcelGDF = gpd.read_file(file).to_crs(merged.crs)

test = get_parcel_area(parcelGDF, 'suss_10005')
dic['suss_10005'] = test

## VIRGINIA
# VA had a statewide dataset...long but can grab in one shot
vacountylist = [county.rstrip('_parcels') for county in counties if '51' in county]
vacountylist.sort()

parcelGDF = gpd.read_file(
    join(
        network_drive,
        project_folder,
        'VA/_statewide/VGIN/Virginia_Parcels_10.1.18/Virginia_Parcel_Dataset_2018Q3.gdb'
        ),
    drive = 'FileGDB',
    layer = 'VA_Parcels_projected').to_crs(merged.crs)

parcelGDF['Shape_Area'] = parcelGDF.geometry.area

joined = centroidGDF.sjoin(
    parcelGDF[['geometry', 'Shape_Area', 'FIPS']],
        how = 'inner',
        predicate = 'within').reset_index()
    # add the county FIPS code for backwards traceability
joined = joined.rename({'FIPS':'county'}, axis = 1)
    
    # we have to use this workaround to get the parcel id's corresponding to smallest parce
dissolved = joined.loc[joined.groupby('index')['Shape_Area'].idxmin()]
dic['virginia'] = dissolved

## PA
pacountylist = [county.rstrip('_parcels') for county in counties if '42' in county]
pacountylist.sort()

# list of county folder names
files = listdir(join(network_drive, project_folder, 'PA'))
# return files that are not in the list of those we've dealt with already
filesToProcess = [file for file in files + pacountylist if file not in pacountylist]
filesToProcess

shapefiles = []
for root, dirs, files in walk(join(network_drive, project_folder, 'PA')):
    # if we have reached files
    if len(files) > 0:
        # find files with 'parcel' and '.shp'
        file = [file for file in files if 'arcel' in file and '.shp' in file]
        if len(file) > 0:
            shapefiles.append(join(root, file[0]))

finalshapefiles = {}
for county in filesToProcess:
    for file in shapefiles:
        if county in file:
            finalshapefiles[county]=file

processedFiles = []

for county, file in finalshapefiles.items():
    print(county)
    print(file)
    parcelGDF = gpd.read_file(file).to_crs(merged.crs)
    temp = get_parcel_area(parcelGDF, county)
    dic[county] = temp
    processedFiles.append(county)
    
# we wont find shapefiles that are part of gdb, so do these individually  
[file for file in filesToProcess if file not in processedFiles]  
parcelGDF = gpd.read_file(
    join(
        network_drive,
        project_folder,
        'PA\hunt_42061\StandardDataset\StandardDataset.gdb'),
    drive = 'FileGDB',
    layer = 'TaxParcels_projected').to_crs(merged.crs)

test = get_parcel_area(parcelGDF, 'hunt_24061')

dic['hunt_24061'] = test

parcelGDF = gpd.read_file(
    join(
        network_drive,
        project_folder,
        'PA\luze_42079\LUZ_CO_DATA_SEP_19\LUZ_CO_PARCELS_SEP_19_projected.shp')
    ).to_crs(merged.crs)

test = get_parcel_area(parcelGDF, 'luze_42079')

dic['luze_42079'] = test

## NEW YORK
countylist = [county.rstrip('_parcels') for county in counties if '36' in county]
countylist.sort()

# list of county folder names
files = listdir(join(network_drive, project_folder, 'NY'))
# return files that are not in the list of those we've dealt with already
filesToProcess = [file for file in files + countylist if file not in countylist]
filesToProcess

# ...no new counties in NY

type(dic['alle_24001'])

newData = gpd.GeoDataFrame(pd.concat(dic.values(), ignore_index = True), crs = merged.crs).rename({'index_right':'index_parcel'}, axis = 1)

# make sure nothing wierd happened while concatenating
len(newData) == sum([len(df) for df in dic.values()])

# some parcels are in multiple county files, so we need to aggregate again
# here we will again use the smalles shape area 
dissolved = newData.loc[newData.groupby('index')['Shape_Area'].idxmin()]

final = merged.merge(
    dissolved.drop(['level_0', 'geometry'], axis = 1),
    how = 'left',
    on = 'index')

# check for 1-to-1 and unique index
len(final) == len(merged)

len(final[final['index'].duplicated()]) == 0

## Write our data out for analysis in R!
final.to_file('./data/solar_analysis_data.geojson', driver = 'GeoJSON')
final.to_csv('./data/solar_analysis_data.csv')
