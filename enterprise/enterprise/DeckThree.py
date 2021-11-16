"""
Created on Fri Jul  2 23:43:40 2021
@author: picard

This module contain analysis functions
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.stats.stats import pearsonr
from scipy.stats import sem, t
from sklearn.metrics import mean_squared_error
from scipy import ndimage
import netCDF4 as nc
from matplotlib import pyplot as plt
################################################################################
class interpolate:
    def regrid(data, out_x, out_y):
        m = max(data.shape[0], data.shape[1])
        y = np.linspace(0, 1.0/m, data.shape[0])
        x = np.linspace(0, 1.0/m, data.shape[1])
        interpolating_function = RegularGridInterpolator((y, x), data)
        yv, xv = np.meshgrid(np.linspace(0, 1.0/m, out_y), np.linspace(0, 1.0/m, out_x))
        return interpolating_function((xv, yv))

class statistics:
    def corr(obs,exp):
        exp=exp.flatten()
        obs=obs.flatten()
        for i in range(len(obs)):
            if np.isnan(obs[i]) == True:
                exp[i]=np.nan
            elif np.isnan(exp[i]) == True:
                obs[i]=np.nan
        exp=exp[~(np.isnan(exp))]
        obs=obs[~(np.isnan(obs))]
        cor=pearsonr(obs, exp)
        return cor[0],cor[1]
    
    def rmse(obs,exp):
        exp=exp.flatten()
        obs=obs.flatten()
        for i in range(len(obs)):
            if np.isnan(obs[i]) == True:
                exp[i]=np.nan
            elif np.isnan(exp[i]) == True:
                obs[i]=np.nan
        exp=exp[~(np.isnan(exp))]
        obs=obs[~(np.isnan(obs))]
        return np.sqrt(mean_squared_error(obs,exp))

    def bias(obs,exp):
        exp=exp.flatten()
        obs=obs.flatten()
        for i in range(len(obs)):
            if np.isnan(obs[i]) == True:
                exp[i]=np.nan
            elif np.isnan(exp[i]) == True:
                obs[i]=np.nan
        exp=exp[~(np.isnan(exp))]
        obs=obs[~(np.isnan(obs))]
        return np.mean(exp-obs)
    
    def con_int(data, confidence = 0.95):
        data=data[~(np.isnan(data))]
        n=len(data)
        std_err=sem(data)
        h=std_err*t.ppf((1+confidence)/2,n-1)
        return(h)
    
    def remove_outlier(data):
        ''' following http://mathcenter.oxford.emory.edu/site/math117/shapeCenterAndSpread/ '''
        data=np.array(data)
        data=data[~(np.isnan(data))]
        q3=np.percentile(data, 75)
        q1=np.percentile(data, 25)
        iqr=q3-q1
        mn=q1-1.5*iqr
        mx=q3+1.5*iqr
        data1=data[((data<=mx)&(data>=mn))]
        return data1

class cyclone:
    def vel2polvel(u,v,lat,lon,lat0=None,lon0=None):
        if len(np.shape(lat))==1:
            lon,lat=np.meshgrid(lon,lat)
        if lat0 is None :
            s=np.sqrt(u**2+v**2)
            lat0=lat[np.where(s==np.nanmin(s))[0][0],np.where(s==np.nanmin(s))[1][0]]
            lon0=lon[np.where(s==np.nanmin(s))[0][0],np.where(s==np.nanmin(s))[1][0]]
        ref_lon=lon-lon0
        ref_lat=lat-lat0
        v_r=np.ones(np.shape(u))*np.nan
        v_th=np.ones(np.shape(u))*np.nan
        for i in range(np.shape(u)[0]):
            for j in range(np.shape(u)[1]):
                v_r[i,j]=(u[i,j]*ref_lon[i,j]+v[i,j]*ref_lat[i,j])/np.sqrt(ref_lon[i,j]**2+ref_lat[i,j]**2)
                #v_theta
                pos=np.sqrt(ref_lon[i,j]**2+ref_lat[i,j]**2)
                th=np.abs(np.arctan(ref_lat[i,j]/ref_lon[i,j]))
                if np.cos(th)<0:
                    cos=0
                else:
                    cos=np.cos(th)
                x_proj=pos/cos
                x=ref_lon[i,j]
                y=ref_lat[i,j]
                #1st quadrent
                if x > 0 and y > 0:
                    v_th[i,j]=(u[i,j]*(x-x_proj)+v[i,j]*y)/np.sqrt((x-x_proj)**2+y**2)
                #2nd quadrent
                if x < 0 and y > 0:
                    v_th[i,j]=(u[i,j]*(-x-x_proj)-v[i,j]*y)/np.sqrt((-x-x_proj)**2+y**2)
                #3rd quadrent
                if x < 0 and y < 0:
                    v_th[i,j]=(u[i,j]*(x+x_proj)+v[i,j]*y)/np.sqrt((x+x_proj)**2+y**2)
                #4th quadrent
                if x > 0 and y < 0:
                    v_th[i,j]=(u[i,j]*(-x+x_proj)-v[i,j]*y)/np.sqrt((-x+x_proj)**2+y**2)
                #Y-axis
                if x==0 :
                    if y>0 :
                        v_th[i,j]=-u[i,j]
                    elif y<0 :
                        v_th[i,j]=u[i,j]
                elif y==0 :
                    if x>0 :
                        v_th[i,j]=v[i,j]
                    elif x<0 :
                        v_th[i,j]=-v[i,j]
        return (v_r,v_th)
    
    def agm_avg(datapath,t=0):
        data=nc.Dataset(datapath)
        pressure_level=data['lev'][:]
        lat=data['lat'][:]
        lon=data['lon'][:]
        slp=data['slp'][t,:,:]
        lon0=lon[np.where(slp==np.nanmin(slp))[1]][0]
        lat0=lat[np.where(slp==np.nanmin(slp))[0]][0]
        lon,lat=np.meshgrid(lon,lat)
        r=np.sqrt((lon-lon0)**2+(lat-lat0)**2)
        r_uniq=np.unique(r)
        r_lim=r_uniq[r_uniq<=2]
        azm_v_r=[]
        azm_v_th=[]
        for k,pres in enumerate(pressure_level):
            u=data['u'][t,k,:,:]
            v=data['v'][t,k,:,:]
            v_r,v_th=cyclone.vel2polvel(u, v, lat, lon, lat0=lat0, lon0=lon0)
            avg_v_th=np.ones(len(r_lim))*np.nan
            avg_v_r=np.ones(len(r_lim))*np.nan
            for i,r_val in enumerate(r_lim):
                avg_v_th[i]=np.nanmean(v_th[r==r_val])
                avg_v_r[i]=np.nanmean(v_r[r==r_val])
            azm_v_r.append(avg_v_r)
            azm_v_th.append(avg_v_th)
        azm_v_r=np.array(azm_v_r)
        azm_v_th=np.array(azm_v_th)
        fig,ax=plt.subplots()
        fig1=ax.contourf(r_lim,pressure_level,azm_v_th,levels=20,cmap='jet')
        cbar=fig.colorbar(fig1)
        sigma_y = 0
        sigma_x = 100
        sigma = [sigma_y, sigma_x]
        intrp_v_r=azm_v_r
        intrp_v_r[np.isnan(intrp_v_r)]=0
        y = ndimage.filters.gaussian_filter(intrp_v_r, sigma, mode='constant')
        cs=ax.contour(r_lim,pressure_level,y,linewidths=0.3,colors='black',levels=range(-40,50,5))
        ax.clabel(cs,cs.levels,inline=True,fontsize=5,fmt='%1.0f')
        print('contour done')
        ax.invert_yaxis()
        ax.set_yscale('log')
        ax.set_xlim(0,1.95)
        ax.set_ylim(1000,100)
        ax.set_yticks(list(range(1000,400,-200))+list(range(500,0,-100)))
        leb=[str(pres)+' hPa' for pres in list(range(1000,400,-200))+list(range(500,0,-100))]
        ax.set_yticklabels(leb)
        ax.set_xlabel('Radius (degree)')
        plt.savefig(str(t)+'.png',dpi=300)
        plt.show()
        plt.clf()