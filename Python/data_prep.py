# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import rasterio as rio
import geopandas as gpd
from os.path import join
import fiona 

network_drive = 'S:\\'
project_folder = 'SolarAnalysis\Data\Vector'

cpk_solar_file = join(network_drive, project_folder, 'Vector\CPK_solarJun21_final.shp')
cpk_states_file = join(network_drive, project_folder, 'Vector\CPK_states.geojson')
soil_file = join(network_drive, project_folder, 'primesoilscbw')

solarGDF = gpd.read_file(cpk_solar_file)
statesGDF = gpd.read_file(cpk_states_file, driver = 'GeoJSON')
statesGDF.head()
statesGDF.to_file(join(network_drive, project_folder, 'Vector\CPK_states.shp'))

parcel_file = join(network_drive, project_folder, 'Vector/parcels_CBW_bycounty.gdb')
# parcel_layers = fiona.listlayers(parcel_file)


    
