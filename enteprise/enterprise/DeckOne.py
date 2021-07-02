"""
Created on Thu Jul  1 15:44:47 2021
@author: picard

This module contains plotting related fuctions

"""

import pandas as pd
import numpy as np
import datetime as dt
import geopandas as gpd
from py3grads import Grads
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.stats.stats import pearsonr
from sklearn.metrics import mean_squared_error
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams["axes.linewidth"] = 3
plt.rcParams["patch.linewidth"] =3

############################ files ###########################################

shape_file_path='/home/picard/maps-master/States/Admin2.shp'

############################## global functions ###############################

def find_res(x1,x2):
    x1=round(x1)
    x2=round(x2)
    for i in [100,50,20,10,5,2,1]:
        if (int(x2-x1))%i==1 and int(x2-x1)==i :
            pass
        elif (int(x2-x1))%i==1:
            break
        elif x2-x1 > i:
            break
        return i
    
#############################################################################

class SpatialPlots:
    
    def initialize(lat1,lat2,lon1,lon2,res=0,lb_size=14,rotate=0,linewidth=0.8,dec_lim=0):
        ax=plt.axes(projection = ccrs.PlateCarree())
        if ((lat1>0 and lat2<45) and (lon1>50 and lon2<110)):
            if lon2-lon1>=5 or lat2-lat1>=5:
                states=gpd.read_file(shape_file_path)
                ax.add_geometries(states.geometry, crs = ccrs.PlateCarree(),facecolor='none', edgecolor='k',linewidth=linewidth) 
        if res==0:
            y_res=find_res(lat1, lat2)
            x_res=find_res(lon1, lon2)
        else:
            y_res=x_res=res
        ax.set_xticks(np.arange(round(lon1,dec_lim),round(lon2,dec_lim)+1,x_res), crs=ccrs.PlateCarree())
        ax.set_yticks(np.arange(round(lat1,dec_lim),round(lat2,dec_lim)+1,y_res), crs=ccrs.PlateCarree())
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        ax.tick_params(axis='x', labelsize=lb_size, rotation=rotate)
        ax.tick_params(axis='y', labelsize=lb_size)
        ax.set_xlim(lon1,lon2)
        ax.set_ylim(lat1,lat2)
        return ax
    
    def initialize_subplots(lat1,lat2,lon1,lon2,plot_no=1,res=0,lb_size=14,rotate=0,linewidth=0.8,dec_lim=0):
        k=10
        if lat2-lat1 > lon2-lon1 :
            key='horizontal'
        else:
            key='vertical'
        sqrs=[1,4,9,16,32]
        one_lines=[2,3]
        two_lines=[6,8,10]
        three_lines=[12]
        if plot_no in sqrs:
            fig_sg=(k*(lon2-lon1)/(lat2-lat1),k*(lat2-lat1)/(lon2-lon1))
            x_nos=y_nos=np.sqrt(plot_no)
        elif plot_no in one_lines:
            if key == 'horizontal' :
                fig_sg=(plot_no*k*(lon2-lon1)/(lat2-lat1),k*(lat2-lat1)/(lon2-lon1))
                x_nos=plot_no
                y_nos=1
            else:
                fig_sg=(k*(lon2-lon1)/(lat2-lat1),plot_no*k*(lat2-lat1)/(lon2-lon1))
                x_nos=1
                y_nos=plot_no
        elif plot_no in two_lines:
            if key == 'horizontal' :
                fig_sg=(int(0.5*plot_no*k*(lon2-lon1)/(lat2-lat1)),2*k*(lat2-lat1)/(lon2-lon1))
                x_nos=plot_no*0.5
                y_nos=2
            else:
                fig_sg=(2*k*(lon2-lon1)/(lat2-lat1)),int(0.5*plot_no*k*(lat2-lat1)/(lon2-lon1))
                x_nos=2
                y_nos=plot_no*0.5
        elif plot_no in three_lines:
            if key == 'horizontal' :
                fig_sg=(int((plot_no/3)*k*(lon2-lon1)/(lat2-lat1)),3*k*(lat2-lat1)/(lon2-lon1))
                x_nos=plot_no/3
                y_nos=3
            else:
                fig_sg=(3*k*(lon2-lon1)/(lat2-lat1)),int((plot_no/3)*k*(lat2-lat1)/(lon2-lon1))
                x_nos=3
                y_nos=plot_no/3
        if res==0:
            y_res=find_res(lat1, lat2)
            x_res=find_res(lon1, lon2)
        else:
            y_res=x_res=res
        ax=[]
        fig=plt.figure(figsize=fig_sg)
        for i in range(plot_no):
            ax.append(fig.add_subplot(int(x_nos),int(y_nos),i+1,projection=ccrs.PlateCarree()))
            if ((lat1>0 and lat2<45) and (lon1>50 and lon2<110)):
                if lon2-lon1>=5 or lat2-lat1>=5:
                    states=gpd.read_file(shape_file_path)
                    ax[i].add_geometries(states.geometry, crs = ccrs.PlateCarree(),facecolor='none', edgecolor='k',linewidth=linewidth) 
            ax[i].set_xticks(np.arange(round(lon1,dec_lim),round(lon2,dec_lim)+1,x_res), crs=ccrs.PlateCarree())
            ax[i].set_yticks(np.arange(round(lat1,dec_lim),round(lat2,dec_lim)+1,y_res), crs=ccrs.PlateCarree())
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax[i].xaxis.set_major_formatter(lon_formatter)
            ax[i].yaxis.set_major_formatter(lat_formatter)
            ax[i].set_xlim(lon1,lon2)
            ax[i].set_ylim(lat1,lat2)
            ax[i].tick_params(axis='x', labelsize=lb_size, rotation=rotate)
            ax[i].tick_params(axis='y', labelsize=lb_size)
        plt.subplots_adjust(left=0.05,
                            bottom=0.05,
                            right=0.95,
                            top=0.95,
                            wspace=0.2,
                            hspace=0.2,)
        return fig,ax    