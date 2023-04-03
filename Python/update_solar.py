# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 13:30:53 2022

@author: mevans
"""

import geopandas as gpd
import pandas as pd
from os.path import join

from shapely.geometry import Point, Polygon
import numpy as np
from numpy.random import normal
import json

ROOT = 'S:\SolarAnalysis\Data\Vector'
fisrt_file = join(ROOT, 'CPK_solarJun21_firstyear.shp')
annual_file = join(ROOT, 'CPK_solarJun21_annual.shp')

# read in our existing solar data
annual = gpd.read_file(annual_file)
first = gpd.read_file(fisrt_file)

# we have new observations for a single year
new_files = [f'C:/Users/mevans/Downloads/solar_{state}2022_Jun21.geojson' for state in ['DE', 'MD', 'NY', 'PA', 'VA']]

gdfs = [gpd.read_file(f, driver = 'GeoJON') for f in new_files]
gdf = gpd.GeoDataFrame(pd.concat(gdfs, ignore_index = True), crs = gdfs[0].crs)
gdf['system_tim'] = '2022-11-01'
gdf['year_right'] = 2022

# update our first year file
new_first = pd.concat([first, gdf.drop(columns = 'system:time_start')],
                      axis = 0,
                      ignore_index = True)

# update our annual file
# step 1: create a complete inventory of arrays present in newest year
last_year = gdf['year_right'].unique() - 1
previous = annual[annual['system_tim'].str.contains('2021')]

# merge our new additions with the existing arrays
updated = pd.concat([previous, gdf[list(set(gdf.columns) & set(annual.columns))]], # only grab the columns gdf has in common with annual
                    axis = 0,
                    ignore_index = True)

new_annual = pd.concat([annual, updated],
                       axis = 0,
                       ignore_index = True)

new_annual[['geometry', 'system_tim']][new_annual['system_tim'].str.contains('2021')].to_file('CPKsolar_2021.geojson')
new_annual[['geometry', 'system_tim']][new_annual['system_tim'].str.contains('2022')].to_file('CPKsolar_2022.geojson')
