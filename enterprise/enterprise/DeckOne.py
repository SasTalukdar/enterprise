"""
Created on Thu Jul  1 15:44:47 2021
@author: picard

This module contains plotting related fuctions

"""
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['font.size']=16
plt.rcParams['xtick.labelsize']=16
plt.rcParams['ytick.labelsize']=16
plt.rcParams["axes.linewidth"] = 3
plt.rcParams["patch.linewidth"] =3


############################ files ###########################################

default_shape_file_path='/home/picard/maps-master/States/Admin2.shp'

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
    
######################## global variables #####################################

labels=['(a)','(b)','(c)','(d)','(e)','(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']
    
#############################################################################

class SpatialPlots:
    def initialize(lat1,lat2,lon1,lon2,res=0,lb_size=14,rotate=0,linewidth=0.8,dec_lim=0,shape_file_path=None,title=None):
        ax=plt.axes(projection = ccrs.PlateCarree())
        if shape_file_path is not None:
            if type(shape_file_path) == str:
                states=gpd.read_file(shape_file_path)
            else:
                states=shape_file_path
            ax.add_geometries(states.geometry, crs = ccrs.PlateCarree(),facecolor='none', edgecolor='k',linewidth=linewidth)
        elif ((lat1>0 and lat2<45) and (lon1>50 and lon2<110)) :
            if (lon2-lon1>=5 or lat2-lat1>=5):
                states=gpd.read_file(default_shape_file_path)
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
        if title is not None:
            ax.set_title(title)
        return ax
    
    def initialize_subplots(lat1,lat2,lon1,lon2,plot_no=1,res=0,lb_size=14,rotate=0,
                            linewidth=0.8,dec_lim=0,shape_file_path=None,add_label=False,titles=None,
                            fig_sg=None,x_nos=None,y_nos=None, coastlines=False):
        plt.rcParams['xtick.major.size']=5
        plt.rcParams['ytick.major.size']=5
        plt.rcParams['xtick.major.width']=3
        plt.rcParams['ytick.major.width']=3
        k=10
        if lat2-lat1 > lon2-lon1 :
            key='vertical'
        else:
            key='horizontal'
        sqrs=[1,4,9,16,32]
        one_lines=[2,3]
        two_lines=[6,8,10]
        three_lines=[12]
        if x_nos is None and y_nos is None:
            if plot_no in sqrs:
                if fig_sg is None:
                    fig_sg=(k*(lon2-lon1)/(lat2-lat1),k*(lat2-lat1)/(lon2-lon1))      # Remove the if condition after algorithm is properly tuned
                x_nos=y_nos=np.sqrt(plot_no)
            elif plot_no in one_lines:
                if key == 'horizontal' :
                    if fig_sg is None:
                        fig_sg=(plot_no*k*(lon2-lon1)/(lat2-lat1),k*(lat2-lat1)/(lon2-lon1))
                    x_nos=plot_no
                    y_nos=1
                else:
                    if fig_sg is None:
                        fig_sg=(k*(lon2-lon1)/(lat2-lat1),plot_no*k*(lat2-lat1)/(lon2-lon1))
                    x_nos=1
                    y_nos=plot_no
            elif plot_no in two_lines:
                if key == 'horizontal' :
                    if fig_sg is None:
                        fig_sg=(int(0.5*plot_no*k*(lon2-lon1)/(lat2-lat1)),2*k*(lat2-lat1)/(lon2-lon1))
                    x_nos=plot_no*0.5
                    y_nos=2
                else:
                    if fig_sg is None:
                        fig_sg=(2*k*(lon2-lon1)/(lat2-lat1)),int(0.5*plot_no*k*(lat2-lat1)/(lon2-lon1))
                    x_nos=2
                    y_nos=plot_no*0.5
            elif plot_no in three_lines:
                if key == 'horizontal' :
                    if fig_sg is None:
                        fig_sg=(int((plot_no/3)*k*(lon2-lon1)/(lat2-lat1)),3*k*(lat2-lat1)/(lon2-lon1))
                    x_nos=plot_no/3
                    y_nos=3
                else:
                    if fig_sg is None:
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
            if shape_file_path is not None:
                if type(shape_file_path) == str:
                    states=gpd.read_file(shape_file_path)
                else:
                    states=shape_file_path
                ax[i].add_geometries(states.geometry, crs = ccrs.PlateCarree(),facecolor='none', edgecolor='k',linewidth=linewidth)
            if ((lat1>0 and lat2<45) and (lon1>50 and lon2<110)):
                if lon2-lon1>=5 or lat2-lat1>=5 and shape_file_path == None:
                    states=gpd.read_file(default_shape_file_path)
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
            if add_label is True:
                ax[i].text(lon1,lat2+(lat2-lat1)/17,labels[i])
            if titles is not None:
                ax[i].set_title(titles[i])
            if coastlines is True:
                ax[i].coastlines()
        plt.subplots_adjust(left=0.05,
                            bottom=0.05,
                            right=0.95,
                            top=0.95,
                            wspace=0.2,
                            hspace=0.2,)
        return fig,ax
    
    def polar_plot(lon,lat,val,lon0=None,lat0=None,r_lim=None,t_lim=None,th_off=np.pi/2.0,t_tick=[0,np.pi*0.5,np.pi,np.pi*1.5],r_tick=None,r_pos=None,c_map='jet'):
        if lon0 is None:
            if len(np.shape(lon))==2:
                lon0=np.nanmean(lon,axis=1)[0]
            else:
                lon0=np.nanmean(lon)
        if lat0 is None:
            if len(np.shape(lat))==2:
                lat0=np.nanmean(lat,axis=0)[0]
            else:
                lat0=np.nanmean(lat)
        r=np.sqrt((lon-lon0)**2+(lat-lat0)**2)
        t=np.arctan2(lon-lon0,lat-lat0)
        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        if r_lim is not None:
            ax.set_ylim(0,r_lim)
        if t_lim is not None:
            ax.set_xlim(0,t_lim)
        ax.contourf(t,r,val,cmap=c_map)
        ax.set_theta_direction(-1)
        # ax.set_theta_zero_location("N")
        ax.set_theta_offset(th_off)
        ax.set_xticks(t_tick)
        if r_tick is not None:
            ax.set_yticks(r_tick)
        if r_pos is not None:
            ax.set_rlabel_position(r_pos)
        return fig, ax
        

class AnalyticalPlots:
    def ContouredScatter(x,y,bins=100,level=30,cmap='jet',savename='contoured_scatter',limits=None):
        xx=np.linspace(min(x),max(x),bins)
        yy=np.linspace(min(y),max(y),bins)
        zz=np.zeros((bins,bins))
        for i in range(len(x)):
            dis=np.abs(xx-x[i])
            ii=np.where(dis==np.nanmin(dis))[0]
            dis=np.abs(yy-y[i])
            jj=np.where(dis==np.nanmin(dis))[0]
            zz[ii,jj]=zz[ii,jj]+1
        plt.contourf(xx,yy,zz,cmap=cmap,level=level)
        if limits is not None:
            plt.xlim(limits[0],limits[1])
            plt.ylim(limits[2],limits[3])
        plt.colorbar()
        plt.savefig(savename+'.png',dpi=300)
        return None
    
    ''' Work under progress 
    def gridded_correl_plots(obs_data,exp_data,k=10):
        if type(obs_data) != list:
            obs_data=[obs_data]
        if type(exp_data) != list:
            exp_data=[exp_data]
        from enterprise.DeckThree import statistics
        ax=[]
        lines=[[1,4,9,16,32],[2,3],[6,8,10],[12]]
        for i in range(len(lines)):
            if i == 0:
                k1=k2=int(np.sqrt(lines[i]))
            elif len(obs_data) in lines[i]:
                k1=int(len(obs_data)/i)
                k2=int(i)
        fig=plt.figure(figsize=(k1*k,k2*k))
        for i in len(obs_data):
            obs=obs_data[i]
            exp=exp_data[i]
            if np.shape(obs) != np.shape(exp):
                from enterprise.DeckThree import interpolate
                exp=interpolate.regrid(exp, np.shape(obs)[0], np.shape(obs)[1])
            rmse=statistics.rmse(obs, exp)
            corr=statistics.corr(obs, exp)
            bias=statistics.bias(obs, exp)
            ax.append(fig.add_subplot(k2,k1,i+1))
            for i in range(3):
                plt.scatter(obs,all_exp[i][dt],alpha=0.4,s=3,color=colors[i])
                z=np.polyfit(obs[dt],all_exp[i][dt],1)
                p=np.poly1d(z)
                if dt == 1:
                    plt.plot(obs[dt],p(obs[dt]),'--',color=colors[i],label=dt_name[i])
                else:
                    plt.plot(obs[dt],p(obs[dt]),'--',color=colors[i])
            ax[dt].set_ylim(0,140/2.3)
            ax[dt].set_xlim(0,140/2.3)
            # ax[dt].set_ylim(0,140)
            # ax[dt].set_xlim(0,140)
            # ax[dt].text(90, 130, 'BIAS  RMSE  CORR ',size=17,c='k')
            # ax[dt].text(90, 123, str(round(all_bias[0][dt],2))+'  '+str(round(all_rmse[0][dt],2))+'  '+str(round(all_cor[0][dt],2)),size=17,c=colors[0])
            # ax[dt].text(90, 116, str(round(all_bias[1][dt],2))+'  '+str(round(all_rmse[1][dt],2))+'  '+str(round(all_cor[1][dt],2)),size=17,c=colors[1])
            # ax[dt].text(90, 109, str(round(all_bias[2][dt],2))+'  '+str(round(all_rmse[2][dt],2))+'  '+str(round(all_cor[2][dt],2)),size=17,c=colors[2])
            # ax[dt].text(10, 143,titles[dt],size=20)
            ax[dt].text(90/2.3, 130/2.3, 'BIAS  RMSE  CORR ',size=17,c='k')
            ax[dt].text(90/2.3, 123/2.3, str(round(all_bias[0][dt],2))+'  '+str(round(all_rmse[0][dt],2))+'  '+str(round(all_cor[0][dt],2)),size=17,c=colors[0])
            ax[dt].text(90/2.3, 116/2.3, str(round(all_bias[1][dt],2))+'  '+str(round(all_rmse[1][dt],2))+'  '+str(round(all_cor[1][dt],2)),size=17,c=colors[1])
            ax[dt].text(90/2.3, 109/2.3, str(round(all_bias[2][dt],2))+'  '+str(round(all_rmse[2][dt],2))+'  '+str(round(all_cor[2][dt],2)),size=17,c=colors[2])
            ax[dt].text(2, 62.3,titles[dt],size=20)
            plt.xlabel('TRMM (mm $day^{-1}$)', size=20)
            if dt == 0:
                plt.ylabel('Model (mm $day^{-1}$)', size=20)
            plt.legend(loc=4)
            plt.savefig('plots/czil_final_'+year+'.png',dpi=300)  '''