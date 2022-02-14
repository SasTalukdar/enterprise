"""
Created on Mon Jul  5 17:49:20 2021
@author: picard

This module contains GIS related functions
"""

import numpy as np
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
from enterprise import HoloDeck

############################ files ###########################################

default_shape_file_path=HoloDeck.__file__[:-12]+'/data/shape_files/Admin2.shp'

###############################################################################

class GeoDataModifications:
    def transform_earth_explorer(data_path,epsg='32645'):
        from subprocess import run,PIPE
        if data_path.endswith('.tif') or data_path.endswith('.TIF'):
            run(['gdalwarp','-s_srs','EPSG:'+epsg,'-t_srs','EPSG:4326',data_path,'modified_'+data_path])
        else:
            ls = run(['ls',data_path], stdout=PIPE).stdout.splitlines()
            if b'projection_modified' not in ls:
                run(['mkdir',data_path+'/projection_modified'])
            for i in range(len(ls)):
                file=str(ls[i])[2:-1]
                if file.endswith('.tif') or file.endswith('.TIF'):
                    run(['gdalwarp','-s_srs','EPSG:'+epsg,'-t_srs','EPSG:4326',data_path+'/'+file,data_path+'/projection_modified/'+file])
        return None
    
class GeoAnalytics:
    def extract_coord(ds):
        nc=ds.RasterXSize
        nr=ds.RasterYSize
        geotransform= ds.GetGeoTransform()
        XOrigin=geotransform[0]
        YOrigin=geotransform[3]
        pixelWidth=geotransform[1]
        pixelHeight=geotransform[5]
        lons=XOrigin+np.arange(0,nc)*pixelWidth
        lats=YOrigin+np.arange(0,nr)*pixelHeight
        return lons,lats
    
    def is_inside(latitude, longitude, shapefile=default_shape_file_path):
        ind = shapefile
        if type(ind) == str:
            ind_geom = unary_union(list(shpreader.Reader(ind).geometries()))
        else:
            ind_geom = unary_union(ind.geometry)
        ind = prep(ind_geom)
        return ind.contains(sgeom.Point(longitude, latitude))
