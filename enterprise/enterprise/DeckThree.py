"""
Created on Fri Jul  2 23:43:40 2021
@author: picard

This module contains analysis functions
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.stats.stats import pearsonr
from scipy.stats import sem, t
from sklearn.metrics import mean_squared_error
from scipy import ndimage
import netCDF4 as nc
from matplotlib import pyplot as plt
import pandas as pd
import xarray as xr
from scipy.ndimage import rotate, shift
from concurrent.futures import ProcessPoolExecutor
################################################################################
class interpolate:
    def regrid(data, out_x, out_y):
        m = max(data.shape[0], data.shape[1])
        y = np.linspace(0, 1.0/m, data.shape[0])
        x = np.linspace(0, 1.0/m, data.shape[1])
        interpolating_function = RegularGridInterpolator((y, x), data)
        yv, xv = np.meshgrid(np.linspace(0, 1.0/m, out_y), np.linspace(0, 1.0/m, out_x))
        return interpolating_function((xv, yv))
    
    def knn_interpolator(data, out_x, out_y):
        x=data.shape[1]
        y=data.shape[0]
        if out_x>x or out_y>y:
            print('output size not valid')
            return None
        out_data=np.ones((out_y,out_x))*np.nan
        x_dif=round(x*0.5/out_x)
        y_dif=round(y*0.5/out_y)
        for i in range(out_y):
            for j in range(out_x):
                i_start=(2*i-1)*x_dif
                j_start=(2*j-1)*y_dif
                if i_start<0:
                    i_start=0
                if j_start<0:
                    j_start=0
                temp=data[i_start:(2*i+1)*x_dif+1,j_start:(2*j+1)*y_dif+1]
                try:
                    out_data[i,j]=np.bincount(temp.flatten()).argmax()
                except:
                    pass
        return out_data

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
    
    def con_int(data, confidence = 0.95, axis=0):
        data=np.array(data)
        if len(np.shape(data)) == 2:
            h=np.zeros(np.shape(data)[1-axis])
            for i in range(np.shape(data)[1-axis]):
                if axis==0:
                    tmp=data[:,i]
                else:
                    tmp=data[i,:]
                tmp=tmp[~(np.isnan(tmp))]
                n=len(tmp)
                std_err=sem(tmp)
                h[i]=std_err*t.ppf((1+confidence)/2,n-1)
        else:   
            data=data[~(np.isnan(data))]
            n=len(data)
            std_err=sem(data)
            h=std_err*t.ppf((1+confidence)/2,n-1)
        return(h)
    
    def remove_outlier(data):
        ''' following http://mathcenter.oxford.emory.edu/site/math117/shapeCenterAndSpread/ '''
        data=np.array(data)
        q3=np.percentile(data, 75)
        q1=np.percentile(data, 25)
        iqr=q3-q1
        mn=q1-1.5*iqr
        mx=q3+1.5*iqr
        data[((data>mx)|(data<mn))]=np.nan
        data1=data
        return data1
    def rot(par):
        [th,exp_cra,obs_cra,x_c,y_c,x_extnd,y_extnd,obs_res,mod_res,del_x,best_fit_criteria,minimum_corr_significance] = par
        bst_corr = -100 # initialization
        bst_p = 1
        transformed = np.nan
        r = np.nan
        temp_exp = shift(exp_cra,[x_c,y_c])
        temp_exp[temp_exp<0.01] = 0
        temp_exp = rotate(temp_exp,th,reshape=False)
        temp_exp[temp_exp<0.01] = 0
        temp_exp1 = shift(temp_exp,[-x_c,-y_c])
        temp_exp1[temp_exp1<0.01] = 0
        for i in range(-round(x_extnd*obs_res/mod_res),round(x_extnd*obs_res/mod_res)+round(del_x/mod_res),round(del_x/mod_res)):
            for j in range(-round(y_extnd*obs_res/mod_res),round(y_extnd*obs_res/mod_res)+round(del_x/mod_res),round(del_x/mod_res)):
                temp_exp = shift(temp_exp1,[i,j])
                temp_exp[temp_exp<0.01] = 0
                temp_exp = interpolate.regrid(temp_exp, np.shape(obs_cra)[0], np.shape(obs_cra)[1])
                if best_fit_criteria == 'corr' :
                    corr,p = statistics.corr(obs_cra, temp_exp)
                    if corr > bst_corr and p <= minimum_corr_significance :
                        transformed = temp_exp
                        r = [i*mod_res,j*mod_res]
                        bst_corr = corr
                        bst_p = p
        return [bst_corr,bst_p,transformed,r,th]
        
    def cra(obs,exp,obs_res=0.1,mod_res=0.018,threshold=5,sep=1,max_shift=1,max_theta=45,del_x=0.05,del_th=5,
            pivot = 'cg_exp',best_fit_criteria = 'corr',minimum_corr_significance = 0.05,n_jobs=1) :
        '''
    
        Parameters
        ----------
        obs : DataArray
            observation matrix extracted through xarray.
        exp : DataArray
            model output matrix extracted through xarray.
        threshold : float, optional
            theshold for CRA detection. The default is 5.
        sep : integer, optional
            separation between two grids to be considered in the same CRA. The default is 1.
        max_shift : float, optional
            maximum shift allowed during matching (degree). The default is 1.
        max_theta : float, optional
            maximum rotation allowed during matching (degree). The default is 45.
        del_x : float, optional
            increments in a single shift (degree). The default is 0.05.
        del_th : float, optional
            angular increaments during rotation (degree). The default is 5.
        pivot : string, optional
            pivot point used for rotation.
            options: 'cg_exp' : center of gravity of the model part of the CRA
                     'cg_cra' : center of gravity of the CRA
            The default is 'cg_exp'.
        best_fit_criteria : string, optional
            type of best fit criteria used for matching.
            options: 'corr' : Pearson correlation coefficient
            The default is 'corr'.
        minimum_corr_significance : float, optional
            Minimum p-value of correlation to be considered a match. The default is 0.05.
        n_jobs : integer, optional
            number of processors to be used. The default is 1.
    
        Returns
        -------
        DataFrame containing
        cra_no : integer
            unique identifier for each CRA
        lat : array
            latitude array for the CRA
        lon : array
            longitude array for the CRA
        observation : matrix
            Observation part of the CRA
        model : matrix
            model part of the CRA
        shifted : matrix
            model output after shifting
        transformed : matrix
            model output after shifting and rotating
        r : [i,j]
            lateral shift
        theta : float
            rotation
        corr : float
            optimum correlation
        p_val : float
            p value for the optimum correlation
        Match : Yes/No
            Is match found for the CRA ?
        err_dis : float
            displacement error
        err_rot : float
            rotational error
        err_vol : float
            volumetric error
        err_pat : float
            pattern error
        tot_err : float
            total error
    
        '''
        
        obs_res = np.mean([obs.lat[1]-obs.lat[0],obs.lon[1]-obs.lon[0]])
        mod_res = np.mean([exp.lat[1]-exp.lat[0],exp.lon[1]-exp.lon[0]])
        
        extnd = round(max_shift/obs_res)
        exp_reg = xr.DataArray(interpolate.regrid(exp.values, np.shape(obs)[0], np.shape(obs)[1]), coords=obs.coords)
        cra = np.max(np.array([obs,exp_reg]),axis=0)
        cra[cra<threshold]=0
        # entity finder
        idt = np.ones(np.shape(cra))*np.nan
        k=1
        for i in range(sep,np.shape(cra)[0]-sep):
            for j in range(sep,np.shape(cra)[1]-sep):
                if cra[i,j] != 0 : 
                    if np.isnan(idt[i,j]) == True :
                        if k not in idt.flatten() :
                            idt[i,j] = k
                        else :
                            idt[i,j] = k+1
                            k = k+1
                    for m in range(-sep,sep+1):
                        for n in range(-sep,sep+1):
                            if not (m == 0 and n == 0) :
                                if cra[i-m,j-n] != 0 :
                                    if np.isnan(idt[i-m,j-n]) == True:
                                        idt[i-m,j-n]=idt[i,j]
                                    else :
                                        l = idt[i-m,j-n]
                                        idt[idt==l]=idt[i,j]
        for i, val in enumerate(np.unique(idt)):
            if not np.isnan(val):
                idt[idt==val]=i
        variables = ['cra_no','lat','lon','observation_cra','original_cra','shifted_cra','transformed_cra','r','theta','corr','p_val','Match']
        cra_details = pd.DataFrame(columns=variables) # dataframe to save CRA details
        idts=np.unique(idt) # unique identities
        idts=idts[~(np.isnan(idts))] # drop nan identity
        for val in idts:
            print(val)
            match = 'No' # Indication for whether a match is found
            locs = np.where(idt==val) # locations of the CRA
            ind = idt[min(locs[0]):max(locs[0])+1,min(locs[1]):max(locs[1])+1] #indices for the boundary of CRA
            obs_part = obs[min(locs[0]):max(locs[0])+1,min(locs[1]):max(locs[1])+1].values # extract the observations contained in the boundary
            obs_cord = obs[min(locs[0]):max(locs[0])+1,min(locs[1]):max(locs[1])+1].coords # extract the coordinates of the observations
            obs_part[ind!=val]=0 # remove the values outside the CRA
            # extend the boundery by the smaller between extend and length of the CRA area
            if extnd > np.shape(obs_part)[0] :
                x_extnd = np.shape(obs_part)[0]
            else :
                x_extnd = extnd
            if extnd > np.shape(obs_part)[1] :
                y_extnd = np.shape(obs_part)[1]
            else :
                y_extnd = extnd
            if x_extnd > y_extnd :
                y_extnd = x_extnd
            else :
                x_extnd = y_extnd   
            cra_lat = np.arange(min(obs_cord['lat'])-x_extnd*obs_res,max(obs_cord['lat'])+x_extnd*obs_res+obs_res,obs_res) # reference lattitude for the CRA
            cra_lon = np.arange(min(obs_cord['lon'])-x_extnd*obs_res,max(obs_cord['lon'])+x_extnd*obs_res+obs_res,obs_res) # reference longitude for the CRA
            exp_part = exp_reg[min(locs[0]):max(locs[0])+1,min(locs[1]):max(locs[1])+1].values # extract the model values contained in the boundary
            exp_part[ind!=val]=0 # remove the values outside the CRA
            if len(exp_part[exp_part!=0]) > 0 and len(obs_part[obs_part!=0]) > 0 :
                obs_cra = np.zeros((np.shape(obs_part)[0]+2*x_extnd,np.shape(obs_part)[1]+2*y_extnd)) # place to store observations
                obs_cra[x_extnd:np.shape(obs_part)[0]+x_extnd,y_extnd:np.shape(obs_part)[1]+y_extnd]=obs_part # store the observations
                exp_org = exp[round(min(locs[0])*len(exp.lat)/len(exp_reg.lat)):round(max((locs[0]))*len(exp.lat)/len(exp_reg.lat))+1,
                              round(min(locs[1])*len(exp.lon)/len(exp_reg.lon)):round(max((locs[1]))*len(exp.lon)/len(exp_reg.lon))+1].values
                if np.shape(exp_part)[0] > 1 and np.shape(exp_part)[1] > 1 :
                    # it is not possible to regrid without at least 2 values
                    exp_part_reg = interpolate.regrid(exp_part, np.shape(exp_org)[0], np.shape(exp_org)[1])
                    exp_org[exp_part_reg==0]=0     
                exp_cra = np.zeros((np.shape(exp_org)[0]+round(2*x_extnd*len(exp.lat)/len(exp_reg.lat)),
                                   np.shape(exp_org)[1]+round(2*y_extnd*len(exp.lon)/len(exp_reg.lon))))
                exp_cra[round(x_extnd*len(exp.lat)/len(exp_reg.lat)):np.shape(exp_org)[0]+round(x_extnd*len(exp.lat)/len(exp_reg.lat)),
                        round(y_extnd*len(exp.lon)/len(exp_reg.lon)):np.shape(exp_org)[1]+round(y_extnd*len(exp.lon)/len(exp_reg.lon))]=exp_org
                # Lagrangian search
                if pivot == 'cg_exp' :
                    y,x = np.meshgrid(np.arange(np.shape(exp_cra)[1])+1,np.arange(np.shape(exp_cra)[0])+1)
                    x_c = np.sum(exp_cra*x)/np.sum(exp_cra)
                    y_c = np.sum(exp_cra*y)/np.sum(exp_cra)
                    x_c = round(np.shape(exp_cra)[0]/2 - x_c)
                    y_c = round(np.shape(exp_cra)[1]/2 - y_c)
                elif pivot == 'cg_cra' :
                    temp_obs = interpolate.regrid(obs_cra, np.shape(exp_cra)[0], np.shape(exp_cra)[1])
                    temp_cra = np.max(np.array([exp_cra,temp_obs]),axis=0)
                    y,x = np.meshgrid(np.arange(np.shape(temp_cra)[1])+1,np.arange(np.shape(temp_cra)[0])+1)
                    x_c = np.sum(temp_cra*x)/np.sum(temp_cra)
                    y_c = np.sum(temp_cra*y)/np.sum(temp_cra)
                    x_c = round(np.shape(temp_cra)[0]/2 - x_c)
                    y_c = round(np.shape(temp_cra)[1]/2 - y_c)
                if best_fit_criteria == 'corr' :
                    opt = -100 # variable to store optimum correlation coeff
                    opt_p = 1 # variable to store optimum p value
                transformed = obs_cra*np.nan # matrix to store optimally shifted and rotated exp CRA
                r = [np.nan,np.nan]
                theta = np.nan
                pars = [[th,exp_cra,obs_cra,x_c,y_c,x_extnd,y_extnd,obs_res,mod_res,del_x,best_fit_criteria,minimum_corr_significance] for th in range(-max_theta,max_theta+del_th,del_th)]
                with ProcessPoolExecutor(max_workers=n_jobs) as executor:
                    results = executor.map(statistics.rot, pars)
                    for result in results:
                        [corr1,p1,transformed1,r1,th1] = result
                        if corr1 > opt and p1 <= minimum_corr_significance :
                            opt = corr1
                            opt_p = p1
                            transformed = transformed1
                            r = r1
                            theta = th1
                            match = 'Yes'
                if match == 'Yes':
                    i,j = r # add only shift
                    i = round(i/mod_res)
                    j = round(j/mod_res)
                    shifted_exp = shift(exp_cra,[i,j])
                    shifted_exp[shifted_exp<0.01] = 0
                    shifted_exp = interpolate.regrid(shifted_exp, np.shape(obs_cra)[0], np.shape(obs_cra)[1]) # interpolate model to observation grid
                    exp_cra = interpolate.regrid(exp_cra, np.shape(obs_cra)[0], np.shape(obs_cra)[1]) # interpolate model to observation grid
                    cra_details = pd.DataFrame(np.insert(cra_details.values, len(cra_details.index), values = [val,cra_lat,cra_lon,obs_cra,exp_cra,shifted_exp,transformed,r,theta,opt,opt_p,match], axis=0))
            else :
                cra_details = pd.DataFrame(np.insert(cra_details.values, len(cra_details.index), values = [val,cra_lat,cra_lon,obs_part,exp_part,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,'No'], axis=0))
        cra_details.columns = variables
        cra_details['err_dis'] = np.nan
        cra_details['err_rot'] = np.nan
        cra_details['err_vol'] = np.nan
        cra_details['err_pat'] = np.nan
        cra_details['tot_err'] = np.nan
        for i, row in cra_details.iterrows():
            if row.Match == 'Yes' :
                sf = np.std(row.original_cra)
                sx = np.std(row.observation_cra)
                r = statistics.corr(row.observation_cra,row.original_cra)[0]
                r_s = statistics.corr(row.observation_cra,row.shifted_cra)[0]
                r_opt = statistics.corr(row.observation_cra,row.transformed_cra)[0]
                f_bar = np.mean(row.transformed_cra)
                x_bar = np.mean(row.observation_cra)
                cra_details['err_dis'][i] = 2*sf*sx*(r_s-r)
                cra_details['err_rot'][i] = 2*sf*sx*(r_opt-r_s)
                cra_details['err_vol'][i] = (f_bar-x_bar)**2
                cra_details['err_pat'][i] = 2*sf*sx*(1-r_opt)+(sf-sx)**2
                cra_details['tot_err'][i] = statistics.rmse(row.observation_cra,row.original_cra)**2
        return cra_details

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
