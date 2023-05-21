#!/usr/bin/env python
# coding: utf-8

# In[2]:


import glob
import sys


from netCDF4 import Dataset
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon



import os 
# choose directory where .nc file is located
os.chdir('E:\\Senitinel 1-2_SoilMoisture_FieldData-Updated\\_BrijMohanBohra\\SMOS\\filter_data\\filter_main_grid')

# folder_path where .ncfile is located
folder_path = os.getcwd()
#     extent a list = [(left_lat,lower_long),(left_lat,upper_long),(right_lat,upper_lat),(right_lat,lower_log)]
extent = [( 75.671970,30.358505), ( 75.671970,31.412555), ( 76.761606,31.412555), ( 76.761606, 30.358505)]




def smos_nc_shp_clip(folder_path,extent, extent_name,min_gridx, max_gridx, min_gridy, max_gridy, pixel_width, extract_lat, extract_long, save_folder, csv_file):
#     import library glob,sys,os, netCDF4(Dataset),numpy as np, geopandas as gpd, matplotlib.pylot as plt, shapely.geometry(Point, Polygon)
# crs = string in proj4 format for WGS-84 which is crs for smos file
# extent a list = [(left_lat,lower_long),(left_lat,upper_long),(right_lat,upper_lat),(right_lat,lower_log)]
#  make new directory in folder name clip_file
#     directory = 'clip_file'
#     path = os.path.join(folder_path, directory)
#     os.mkdir(path) 
#     print(os.getcwd())
    # Ignore warning about missing/empty geometries
    import warnings
    warnings.filterwarnings('ignore', 'GeoSeries.notna', UserWarning)
    dir_list = os.listdir(folder_path)
    crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs' 
#     Clip geodataframe for your extent
    clip = gpd.GeoDataFrame()
    clip['geometry'] = None # create a new colum named geometry
    polygon = Polygon(extent) # Create a Shapely polygon from the coordinate-tuple list
    clip.loc[0, 'geometry'] = polygon # Insert the polygon into 'geometry' -column at index 0
    clip.loc[0, 'location'] = extent_name  # Add a new column and insert name of extent
    clip.crs = crs # Set the GeoDataFrame's coordinate system to WGS84 (i.e. epsg code 4326) or using gdf.crs
# generate dataframe for extracting values
    values = pd.DataFrame()
    i = 0
    
    
    for file in dir_list:
#         import netcdf file
        nc_file = Dataset(file)
#     extract variables
        lat = nc_file.variables['Latitude'][:]
        long = nc_file.variables['Longitude'][:]
        sm = nc_file.variables['Soil_Moisture'][:] 
#         make pandadataframe outof variables
        col = {'lats': lat,
        'longs': long,
        'soil_moisture': sm}
        nc_file_df = pd.DataFrame(col)
        nc_file_sub_df = nc_file_df.dropna()
        nc_file_sub_df.reset_index(drop= True, inplace=True)
#         make geodataframe out of dataframe using lat and lon column and crs file 
        gdf = gpd.GeoDataFrame(nc_file_sub_df, geometry=gpd.points_from_xy(nc_file_sub_df.longs, nc_file_sub_df.lats), crs=crs)
#         print(gdf.shape)
# Clip the data using GeoPandas clip
        points_clip = gpd.clip(gdf, clip)
        if points_clip.shape[0] >10 :
#         if not points_clip.empty:
#             for krigging interpoaltion
            clip_lats= points_clip.iloc[:,0]
            clip_long = points_clip.iloc[:,1]
            clip_SM = points_clip.iloc[:,2]
            OK = OrdinaryKriging(
                    clip_long,
                    clip_lats,
                    clip_SM,
                    variogram_model = "linear",
                    verbose =False,
                    enable_plotting = False,
                    coordinates_type= 'geographic'
                    )
            UK = UniversalKriging(
                    clip_long,
                    clip_lats,
                    clip_SM,
                    variogram_model="linear",
                    drift_terms=["regional_linear"])
            gridx = np.arange(min_gridx,max_gridx,0.01)
            gridy = np.arange(min_gridy,max_gridy,0.01)
            z1, ss1 = OK.execute("grid", gridx, gridy)
            z2,ss2 = UK.execute("grid",gridx,gridy)
#             extracting values
            upper_left_long = gridx[0]
            upper_right_lat = gridy[len(gridy)-1]
            pixel_position_x = round((extract_long-upper_left_long)/pixel_width-0.5)
            pixel_position_y = round((upper_right_lat- extract_lat)/pixel_width-0.5)
#   fill extracting values in values dataframe with other attributes
#             i stands for each file number from start.
#             print(i)
#             print(extract_long, extract_lat)
#             print(upper_left_long, upper_right_lat)
#             print(pixel_position_y, pixel_position_x)
            values.loc[i,'date'] = file[19:27]
            values.loc[i,'lat'] = extract_lat
            values.loc[i,'long'] = extract_long
            values.loc[i,'ok_krig'] = z1[pixel_position_y][pixel_position_x]
            values.loc[i,'uk_krig'] = z2[pixel_position_y][pixel_position_x]
            i = i+1
            
#             print(file[19:27],z1[pixel_position_y][pixel_position_x],z2[pixel_position_y][pixel_position_x])
            

#             # Determine the output path for the Shapefile
#             outfp = file[:-2]+"shp"
#     # Write the data into that Shapefile
# #             savefile = os.path.join(path,outfp)
#             points_clip.to_file(outfp)         
        else:
            pass
        
    #do something
    path = os.path.join(save_folder,csv_file)
    values.to_csv(path)
   

        
        
        
        
           
os.chdir('E:\\Senitinel 1-2_SoilMoisture_FieldData-Updated\\_BrijMohanBohra\\SMOS\\filter_data\\filter_buffer_grid')
folder_path = os.getcwd()
extent = [( 75.671970,30.358505), ( 75.671970,31.412555), ( 76.761606,31.412555), ( 76.761606, 30.358505)]
extent_name = "urna_extent"
min_gridx = 75.6
max_gridx = 76.8
min_gridy = 30.3
max_gridy = 31.5
pixel_width = 0.01

save_folder = 'E:\\Senitinel 1-2_SoilMoisture_FieldData-Updated\\_BrijMohanBohra\\SMOS\\filter_data\\extract_csv'
buffer_grid_path = os.path.join(save_folder,'buffer_grid.csv')
df = pd.read_csv(buffer_grid_path)
df.head()
for j in range(df.shape[0]):
    extract_lat = df.loc[j,'Y']
    extract_long = df.loc[j,'X']
    grid_id = df.loc[j,'Grid_No']
    csv_file = grid_id + '.csv'
    print(extract_lat, extract_long, grid_id, csv_file)
    smos_nc_shp_clip(folder_path, extent, extent_name, min_gridx, max_gridx, min_gridy, max_gridy, pixel_width, extract_lat, extract_long,save_folder,csv_file )
    
         
         


