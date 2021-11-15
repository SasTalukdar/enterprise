
"""
Created on Thu Jul  1 15:44:47 2021
@author: picard

This module contains WRF related functions
"""

import numpy as np
from math import prod
from netCDF4 import Dataset
from osgeo import gdal
from enterprise.DeckFour import GeoAnalytics

######################### global functions ####################################

def relhum(t,p,mxr):
    svp=0.622*np.exp(17.67*(t-273.15)/(t-29.65))
    smr=0.622*svp/(p-svp)
    relh=mxr/(smr*100)
    return relh

def dew(t,rh):
    d=t-(100-rh)/5
    return d

def extract_values(line):
    ind=[]
    values=[]
    ind.append(line.find('='))
    for i in range(len(line)):
        if line[i] == ',' :
            ind.append(i)
    for i in range(len(ind)-1):
        values.append(float(line[ind[i]+1:ind[i+1]]))
    return values
######################## global variables #####################################

labels=['(a)','(b)','(c)','(d)','(e)','(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']

##############################################################################

class wps:
    ''' This class contains scripts related to WPS '''
    def namelist_plot(namelist_path):
        from DeckOne.SpatialPlots import initialize
        with open('namelist.wps','r') as file:
            lines=file.readlines()
            file.close()
        for line in lines:
            if 'max_dom' in line:
                eq=line.find('=')
                c=line.find(',')
                dom_no=int(line[eq+1:c])
            if 'parent_grid_ratio' in line:
                pg_rat=extract_values(line)
            if 'i_parent_start' in line:
                i_start=extract_values(line)
            if 'j_parent_start' in line:
                j_start=extract_values(line)
            if 'e_we' in line:
                we=extract_values(line)
            if 'e_sn' in line:
                sn=extract_values(line)
            if 'dx' in line:
                dx=extract_values(line)[0]
            if 'dy' in line:
                dy=extract_values(line)[0]
            if 'ref_lat' in line:
                ref_lat=extract_values(line)[0]
            if 'ref_lon' in line:
                ref_lon=extract_values(line)[0]
        lon1=ref_lon-(dx*we[0]/222222)
        lon2=ref_lon+(dx*we[0]/222222)
        lat1=ref_lat-(dy*sn[0]/222222)
        lat2=ref_lat+(dy*sn[0]/222222)
        ax=initialize(lat1,lat2,lon1,lon2)
        ax.set_title('WPS Grids')
        lons1=[lon1]
        lons2=[lon2]
        lats1=[lat1]
        lats2=[lat2]
        for k in range(1,dom_no):
            lons1.append(lons1[k-1]+i_start[k]*dx/(prod(pg_rat[:k+1])*111111))
            lons2.append(lons1[k]+we[k]*prod(pg_rat[:k+1])/111111)
            lats1.append(lats1[k-1]+j_start[k]*dy/(prod(pg_rat[:k+1])*111111))
            lats2.append(lats1[k]+sn[k]*prod(pg_rat[:k+1])/111111)
            ax.plot([lons1[k],lons2[k],lons2[k],lons1[k],lons1[k]],[lats1[k],lats1[k],lats2[k],lats2[k],lats1[k]])
        return ax
    
    def inject_lulc(input_lulc,exchange_pairs_arr,output_lulc='geo_em.d01.nc'):
        ''' Give the prepared lulc map in tiff format as input.
            Exchange pairs_arr should be an array of to be interchanged values prefererably an
            array of tuples, e.g. () '''
        data=input_lulc
        ds=gdal.Open(data)
        lulc=ds.GetRasterBand(1).ReadAsArray()
        lon,lat=GeoAnalytics.extract_coord(ds)
        lat=np.flip(lat)
        lulc=np.flipud(lulc)       
        for pair in exchange_pairs_arr:
            lulc[lulc==pair[0]] = pair[1]       
        nc = Dataset(output_lulc,'r+')
        lat2=np.squeeze(nc.variables['XLAT_M'])
        lon2=np.squeeze(nc.variables['XLONG_M'])
        lu_val= np.squeeze(nc.variables['LU_INDEX'][:])
        lon_max=np.max(lon2)
        lon_min=np.min(lon2)
        lat_max=np.max(lat2)
        lat_min=np.min(lat2)
        for i in range(len(lon)):
            if lon[i]<=lon_max and lon[i]>=lon_min:
                for j in range(len(lat)):
                    if lat[j]<=lat_max and lat[j]>=lat_min:
                        dis=np.abs(lon2[0,:]-lon[i])
                        ii=np.where(dis==np.nanmin(dis))[0][0]
                        dis=np.abs(lat2[:,0]-lat[j])
                        jj=np.where(dis==np.nanmin(dis))[0][0]
                        lu_val[jj,ii]=lulc[j,i]
        lu_val=np.expand_dims(lu_val, axis=0)
        nc.variables['LU_INDEX'][:] = lu_val
        print('injection complete')
        return lu_val

class assimilation:
    ''' This class contains scripts related to data assimilation '''
    def insat2littler(data_paths):
        ''' enter date_str in the form of yyyymmddhh and enter a list of data files'''
        def form(a):
            for i,x in enumerate(a):
                for j,y in enumerate(x):
                    a[i][j]=format(y, '20.5f')
            return a
        def extract_date(path):
            mon_list=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
            times=[]
            for strng in path:
                times.append(strng[-20:-16]+str(int(mon_list.index(strng[-23:-20])+1))+strng[-25:-23]+strng[-15:-11]+'00')
            return times
        dat=extract_date(data_paths)
        with open('ob.'+dat[0],'w') as f: 
            for ii,dts in enumerate(data_paths):
                data=Dataset(dts)
                lat_rec=np.array(np.nanmean(data["Latitude"][:],axis=1))
                lon_rec=np.array(np.nanmean(data["Longitude"][:],axis=0))
                lat1=data["Latitude"][:]
                lon1=data["Longitude"][:]
                latt1=lat1.tolist()
                lont1=lon1.tolist()
                latn=form(latt1)
                lonn=form(lont1)
                lev=np.asarray(data["plevels"])
                lev=lev.astype(int)      
                pw_pr=np.squeeze(np.asarray(data["L2_PREC_WATER"]))
                pw_pr[pw_pr==-999]=-888888
                
                temp_pr=np.squeeze(np.asarray(data["TAirPhy"]))
                humd=np.squeeze(np.asarray(data["H2OMMRPhy"]))
                
                humd[humd==-999]=np.nan
                temp_pr[temp_pr==-999]=np.nan
                tempc=temp_pr-273

                z_data=[]
                rh=np.zeros((len(lev),len(lat_rec),len(lon_rec)))
                td1=np.zeros((len(lev),len(lat_rec),len(lon_rec)))
                for i in range(len(lat_rec)):
                    for j in range(len(lon_rec)):
                        rh[:,i,j]=relhum(temp_pr[:,i,j],lev,humd[:,i,j])
                        z_data.append(rh) 
                        td1[:,i,j]=(dew(tempc[:,i,j],rh[:,i,j]))
                td=td1+273
                td=td.astype(int)
                rh=rh.astype(int)
                temp_pr=temp_pr.astype(int)
                td[td==-9223372036854775808]=-888888
                temp_pr[temp_pr==-9223372036854775808]=-888888
                rh[rh==-9223372036854775808]=-888888
                ID='%40s' % ('AIRSRET')
                name='%-40s' % ('get data information')
                p_fm='%-40s' % ('FM-133 AIRSRET')
                src='%40s' % ('SOURCE')
                ele="{:20.5f}".format(-888888.00000) ## ele
                vf="{:10d}".format(5) ### valid fields
                err="{:10d}".format(0) ### error
                warn="{:10d}".format(0) ### warning
                sq_no=1
                s_no="{:10d}".format(sq_no)
                dup="{:10d}".format(0)
                sud='{0:>10}'.format('T') ## padded with 4 space
                bog='{0:>10}'.format('F')
                dis='{0:>10}'.format('F')
                u_time="{:10d}".format(-888888)
                j_day="{:10d}".format(-888888)
                date='%20s' % (dat[ii])
                slp=-888888.00000
                ceil=-777777.00000
                qc=0
                end_1c="{:13.5f}".format(ceil)
                end_2c=('%7d%13.5f' % (qc,ceil))
                qcc="{:7d}".format(qc)
                tail_value="{:7d}".format(5)
                tail_er="{:7d}".format(0)
                tail_warn="{:7d}".format(0)
                x = ('%7d%13.5f' % (qc,slp))
                y = ('%13.5f' % (slp))
                for i in range(len(lat_rec)):
                    for j in range(len(lon_rec)):
                        row1 = (latn[i][j],lonn[i][j],ID,name,p_fm,src,ele,vf,err,warn,s_no,dup,sud,bog,dis,u_time,j_day,date,y,x,x,x,x,x,x,x,x,x,x,x,x,qcc)
                        er = (end_1c,end_2c,x,x,x,x,x,x,x,x,qcc)
                        tail = (tail_value,tail_er,tail_warn)
                        rows=[]
                        for k in range(len(lev)):
                            if (temp_pr[k,i,j] == slp and td[k,i,j] == slp) and rh[k,i,j] == slp :
                                pass
                            else:
                                lv=("{:13.5f}".format(lev[k]*100))
                                temp_prr=('%7d%13.5f' % (qc,temp_pr[k,i,j]))
                                rhh=('%7d%13.5f' % (qc,rh[k,i,j]))
                                tdd=('%7d%13.5f' % (qc,td[k,i,j]))
                                row = (lv,x,temp_prr,tdd,x,x,x,x,rhh,x,qcc)
                                sq_no += 1
                                rows.append(row)
                        if len(rows)>0:
                            f.writelines(row1)
                            f.write("\n")
                            for row in rows:
                                f.writelines(row)
                                f.write("\n")
                            f.writelines(er)
                            f.write("\n")
                            f.writelines(tail)
                            f.write("\n")
        return f
        
    def show_assimilated_points(diag_path,namelist_path=None,res=0,s=np.pi/30,alpha=0.7,lb_size=11):
        from enterprise.DeckOne import SpatialPlots
        initialize_subplots=SpatialPlots.initialize_subplots 
        if type(diag_path) == str:
            diag_path=[diag_path]
        if namelist_path != None:
            with open(namelist_path,'r') as file:
                lines=file.readlines()
                file.close()               
            for line in lines:
                if 'IDD' in line:
                    eq=line.find('=')
                    c=line.find(',')
                    dom_no=int(line[eq+1:c])
                if 'MAXNES' in line:
                    pg_rat=extract_values(line)
                if 'NESTI' in line:
                    i_start=extract_values(line)
                if 'NESTJ' in line:
                    j_start=extract_values(line)
                if 'NESTJX' in line:
                    we=extract_values(line)
                if 'NESTIX' in line:
                    sn=extract_values(line)
                if 'DIS' in line:
                    dx=extract_values(line)[0]
                if 'DIS' in line:
                    dy=extract_values(line)[0]
                if 'PHIC' in line:
                    ref_lat=extract_values(line)[0]
                if 'XLONC' in line:
                    ref_lon=extract_values(line)[0]
            lon1=ref_lon-(dx*we[0]/222.222)
            lon2=ref_lon+(dx*we[0]/222.222)
            lat1=ref_lat-(dy*sn[0]/222.222)
            lat2=ref_lat+(dy*sn[0]/222.222)
            fig,ax=initialize_subplots(lat1,lat2,lon1,lon2,plot_no=len(diag_path),res=res)
            lons1=[lon1]
            lons2=[lon2]
            lats1=[lat1]
            lats2=[lat2]
            for axx in ax:
                for k in range(1,dom_no):
                    lons1.append(lons1[k-1]+i_start[k]*dx/(prod(pg_rat[:k+1])*111111))
                    lons2.append(lons1[k]+we[k]*prod(pg_rat[:k+1])/111111)
                    lats1.append(lats1[k-1]+j_start[k]*dy/(prod(pg_rat[:k+1])*111111))
                    lats2.append(lats1[k]+sn[k]*prod(pg_rat[:k+1])/111111)
                    axx.plot([lons1[k],lons2[k],lons2[k],lons1[k],lons1[k]],[lats1[k],lats1[k],lats2[k],lats2[k],lats1[k]])
            for j, event in enumerate(diag_path):
                with open(event) as f:
                    data=f.readlines()
                f.close()
                in_lons=[]
                in_lats=[]                   
                for i, dia in enumerate(data):
                    if 'OUT DOMAIN' not in dia:
                        in_lats.append(float(dia[48:54]))
                        in_lons.append(float(dia[57:63]))
                ax[j].scatter(in_lons,in_lats,s=s,alpha=alpha)
                ax[j].text(lon1+0.1, lat2+0.5, labels[j],size=lb_size)
        else:
            lats=[]
            lons=[]
            lat1=[]
            lat2=[]
            lon1=[]
            lon2=[]
            for j, event in enumerate(diag_path):
                with open(event) as f:
                    data=f.readlines()
                f.close()
                in_lons=[]
                in_lats=[]                   
                for i, dia in enumerate(data):
                    if 'OUT DOMAIN' not in dia:
                        in_lats.append(float(dia[48:54]))
                        in_lons.append(float(dia[57:63]))
                lats.append(in_lats)
                lons.append(in_lons)
                lat1.append(min(in_lats))
                lat2.append(max(in_lats))
                lon1.append(min(in_lons))
                lon2.append(max(in_lons))
            fig,ax=initialize_subplots(min(lat1),max(lat2), min(lon1),max(lon2),plot_no=len(diag_path),res=res)
            for j in range(len(ax)):
                ax[j].scatter(lons[j],lats[j],s=s,alpha=alpha)
                ax[j].text(min(lon1)+0.1, max(lat2)+0.5, labels[j],size=lb_size)
        return ax