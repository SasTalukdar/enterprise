import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

ocean_shp_fname = shpreader.natural_earth(resolution='50m',category='physical', name='ocean')
ocean_geom = unary_union(list(shpreader.Reader(ocean_shp_fname).geometries()))
ocean = prep(ocean_geom)

def is_ocean(longitude, lattitude):
    return ocean.contains(sgeom.Point(longitude, lattitude))

from py3grads import Grads
import numpy as np
import datetime as dt
from pandas import date_range
import pandas as pd
#from tqdm import tqdm
#from tg_tqdm import tg_tqdm

ga=Grads(verbose=False)

def t(m,d,h):
    if m==5:
        t1=0
    elif m==6:
        t1=31
    t=t1*24+(d-1)*24+h+1
    return t
def tgp(m,d,h):
    if m==4:
        t1=0
    elif m==5:
        t1=30
    elif m==6:
        t1=61
    t=t1*24+(d-1)*24+h+1
    return t

def max_wind(lon,lat,y,m,d,h,res):
    ga('reinit')
    ga('sdfopen /home/picard/sas/mov/data/ecmwf/u_v_10m/era5_hourly_srf'+str(y)+'.nc')
    t2=t(m,d,h)
    ga('set t '+str(t2-run)+' '+str(t2+run))
    ga('set lat '+str(lat-4)+' '+str(lat+4))
    ga('set lon '+str(lon-4)+' '+str(lon+4))
    u=ga.exp('u10')
    u=np.nanmean(u,axis=2)
    v=ga.exp('v10')
    v=np.nanmean(v,axis=2)
    s=np.sqrt(u*u+v*v)
    ga('set t '+str(t2))
    la=ga.exp('lat')
    lo=ga.exp('lon')
    max_w=0
    for i in range(len(s)):
        for j in range(len(s)):
            dist=np.sqrt((i-(len(s)-1)/2)**2+(j-(len(s)-1)/2)**2)
            dist=dist*res
            if(dist<=400):
                if s[i][j]>max_w:
                    max_w=s[i][j]
                    exp_lon=lo[i][j]
                    exp_lat=la[i][j]
                    max_dist=dist
    return [max_w,exp_lon,exp_lat,max_dist]


#run=int(input('Enter window for running mean (+-t): '))
run=6
res=27.75
''' Put the resolution of the data in km here '''
year_rec={}
max_w=[]
for year in range(1999, 2000):
    imd_data=pd.read_csv('/home/picard/sas/mov/data/ecmwf/imd_files/'+str(year)+'.csv')
    if len(imd_data)>0:
        print('Calculation running for '+str(year))
        year_rec[year]={}
        rec_cord=[] # Previous record
        recorded=0
        for mon in [5, 6]:
            if (mon==5):
                m=1
                n=32
                month='May'
            elif (mon==6):
                m=1
                n=31
                month='June'
            print(month)
            for day in range(m,n):
                for hour in range(24):
                    ref_date=dt.datetime(year,mon,day,hour)
                    z_all=[]
                    ga('reinit')
                    ga('sdfopen /home/picard/sas/mov/data/ecmwf/gp/era5_hourly_gp_850_'+str(year)+'.nc')
                    ga('set lat 0 35')
                    ga('set lon 40 90')
                    for dates in date_range(start=ref_date-dt.timedelta(hours=504),end=ref_date+dt.timedelta(hours=1)):
                        temp_mon=dates.month
                        temp_day=dates.day
                        temp_hour=dates.hour
                        time=tgp(temp_mon,temp_day,temp_hour)
                        ga('set t '+str(time))
                        z=ga.exp('z')
                        z_all.append(z)
                    z_all=np.array(z_all)
                    z_mean=np.nanmean(z_all,axis=0)
                
                    u=[]
                    v=[]
                    vor=[]
                    z=[]
                    t2=t(mon,day,hour)
                    for time in range(t2-run,t2+run+1):
                        ga('reinit')
                        ga('sdfopen /home/picard/sas/mov/data/ecmwf/u_v_gp/era5_hourly_u_v_gp_850_'+str(year)+'.nc')
                        ga('set t '+str(time))
                        ga('vor=hcurl(u,v)')
                        ga('set lat 0 35')
                        ga('set lon 40 90')
    
                        u.append(ga.exp('u'))
                        v.append(ga.exp('v'))
                        vor.append(ga.exp('vor'))
                        z.append(ga.exp('z'))
                
                    lat=ga.exp('lat')
                    lon=ga.exp('lon')
                    u=np.array(u)
                    u=np.nanmean(u,axis=0)
                    v=np.array(v)
                    v=np.nanmean(v,axis=0)
                    vor=np.array(vor)
                    vor=np.nanmean(vor,axis=0)
                    z=np.array(z)
                    z=np.nanmean(z,axis=0)
                
                    s=np.sqrt(u*u+v*v)
    
                    ''' We don't need negative vorticities. Convert them into nan values '''
                    vor[vor<=0.00001]=np.nan
                
                    ''' Find the coordinates of the local vorticity maximas. '''
                    max_vor=[]
                    for row in range(1,np.shape(vor)[0]-1):
                        for col in range(1,np.shape(vor)[1]-1):
                            x=vor[row][col]
                            if ((x>=vor[row-1][col]) and (x>=vor[row+1][col])):
                                if ((x>=vor[row][col-1]) and (x>=vor[row][col+1])):
                                    if ((x>=vor[row-1][col-1]) and (x>=vor[row+1][col+1])):
                                        if((x>=vor[row-1][col+1]) and (x>=vor[row+1][col-1])):
                                            max_vor.append([row,col])
                    ''' Record the indices instead of lat, lon values. They can be retrived later '''
                
                    ''' max indices for 400 km radius area '''
                    row_max=np.shape(z)[0]-14
                    col_max=np.shape(z)[1]-14
                
                    ''' Finding the eye '''
                    min_wnd=[]
                    for co_ord in max_vor :
                        row=co_ord[0]
                        col=co_ord[1]
                        if ((row>=200/res and row<=row_max) and (col>=200/res and col<=col_max)) :
                            s_near=s[row-int(200/res):row+int(200/res)+1,col-int(200/res):col+int(200/res)+1]
                            pos_s=np.where(s_near==np.nanmin(s_near))
                            min_wnd.append([pos_s[0][0]+row-int(200/res),pos_s[1][0]+col-int(200/res)])
                
                    z_anm=z-z_mean
                
                    ''' Record indices of local maxima vortices with negative slp within 200 km area '''
                    neg_slp=[]
                    for co_ord in min_wnd :
                        row=co_ord[0]
                        col=co_ord[1]
                        marker=0
                        for row_z in range(row-int(200/res),row+int(200/res)+1):
                            for col_z in range(col-int(200/res),col+int(200/res)+1):
                                if (((row>=int(200/res)) and (row<=row_max)) and ((col>=int(200/res)) and (col<=col_max))):
                                    if (z_anm[row_z][col_z]<0):
                                        dist=np.sqrt((row-row_z)**2+(col-col_z)**2)
                                        dist=dist*res
                                        if dist<=200:
                                            marker=1
                        if( marker==1) :
                            neg_slp.append(co_ord)
                
                    ''' Find the points with surface wind speeds greater than 8.5 m/s within 400 km north '''            
                    ga('reinit')
                    ga('sdfopen /home/picard/sas/mov/data/ecmwf/u_v_10m/era5_hourly_srf'+str(year)+'.nc')
                    ga('set lat 0 35')
                    ga('set lon 40 90')
                    ga('set t '+str(t2-run)+' '+str(t2+run))
    
                    u10=ga.exp('u10')
                    u10=np.nanmean(u10,axis=2)
                    v10=ga.exp('v10')
                    v10=np.nanmean(v10,axis=2)
                    s10=np.sqrt(u10*u10+v10*v10)
                
                    srf_wnd=[]
                    for co_ord in neg_slp :
                        row=co_ord[0]
                        col=co_ord[1]
                        marker=0
                        for row_s in range(row-int(400/res),row+int(200/res)+1):
                            ''' Due to LLJ, scanning radius in the southern sector is reduced '''
                            for col_s in range(col-int(400/res)+1,col+int(400/res)):
                                if (((row>=400/res-1) and (row<=row_max)) and ((col>=400/res-1) and (col<=col_max))):
                                    if (s10[row_s][col_s]>8.5):
                                        dist=np.sqrt((row-row_s)**2+(col-col_s)**2)
                                        dist=dist*res
                                        if dist<=400:
                                            marker=1
                        if( marker==1) :
                            srf_wnd.append(co_ord)
    
                    # My logic using wind direction
                    wnd_dir=[]
                    for co_ord in srf_wnd :
                        row=co_ord[0]
                        col=co_ord[1]
    
                        if (((row>=6) and (row<=row_max)) and ((col>=6) and (col<=col_max))):
                            if ((np.sign(u[row-6][col])==1) and (np.sign(u[row+6][col])==-1)) :
                                if ((np.sign(v[row][col-6])==-1) and (np.sign(v[row][col+6])==1)):
                                    if ((np.sign(u[row-4][col-4])==1) and (np.sign(v[row-4][col-4])==-1)):
                                        if ((np.sign(u[row-4][col+4])==1) and (np.sign(v[row-4][col+4])==1)):
                                            if ((np.sign(u[row+4][col-4])==-1) and (np.sign(v[row+4][col-4])==-1)):
                                                if ((np.sign(u[row+4][col+4])==-1) and (np.sign(v[row+4][col+4])==1)):
                                                    wnd_dir.append(co_ord)
    
                        if (((row>=8) and (row<=row_max)) and ((col>=8) and (col<=col_max))):
                            if ((np.sign(u[row-8][col])==1) and (np.sign(u[row+8][col])==-1)) :
                                if ((np.sign(v[row][col-8])==-1) and (np.sign(v[row][col+8])==1)):
                                    if ((np.sign(u[row-6][col-6])==1) and (np.sign(v[row-6][col-6])==-1)):
                                        if ((np.sign(u[row-6][col+6])==1) and (np.sign(v[row-6][col+6])==1)):
                                            if ((np.sign(u[row+6][col-6])==-1) and (np.sign(v[row+6][col-6])==-1)):
                                                if ((np.sign(u[row+6][col+6])==-1) and (np.sign(v[row+6][col+6])==1)):
                                                    wnd_dir.append(co_ord)
                
                        # In the above lines, it has been taken into notice that the direction of row
                        # and the direction of lattitude are opposite.
                
                    # Remove the duplicates 
                    for i in range(len(wnd_dir)-1):
                        for j in range(i+1,len(wnd_dir)):
                            if (wnd_dir[i]!='NaN') and (wnd_dir[j]!='NaN') :
                                dist=np.sqrt((wnd_dir[i][0]-wnd_dir[j][0])**2+(wnd_dir[i][1]-wnd_dir[j][1])**2)
                                dist=res*dist
                                if ((dist<=400) and (s[wnd_dir[i][0],wnd_dir[i][1]]>=s[wnd_dir[j][0],wnd_dir[j][1]])):
                                    wnd_dir[i]='NaN'
                                elif ((dist<=400) and (s[wnd_dir[i][0],wnd_dir[i][1]]<s[wnd_dir[j][0],wnd_dir[j][1]])):
                                    wnd_dir[j]='NaN'
    
                    lats=[]
                    lons=[]
                    for elements in wnd_dir :
                        if (elements!='NaN') :
                            row=elements[0]
                            col=elements[1]
                            latt=lat[row][0]
                            long=lon[0][col]
                            if (((latt>5) and (latt<26)) and ((long>50) and (long<80))) :
                                if is_ocean(long,latt)==True:
                                    lats.append(latt)
                                    lons.append(long)
                    ''' Only the required lat and lon values are recorded '''
                    if len(lats)>0:
                        year_rec[year][dt.datetime(year,mon,day,hour)]=[lons, lats]
                    '''
                    hist_marker=0
                    latss=[]
                    lonss=[]
                    re_cord_lon=[]
                    re_cord_lat=[]
                    for i in range(len(lats)):
                        for j in range(len(rec_cord)):
                            dist=np.sqrt((lats[i]-rec_cord[j][0])**2+(lons[i]-rec_cord[j][1])**2)
                            dist=dist*res
                            if (dist<=36) :
                                hist_marker=1
                                latss.append(lats[i])
                                lonss.append(lons[i])
                                re_cord_lat.append(rec_cord[j][0])
                                re_cord_lon.append(rec_cord[j][1])
                
                    # This condition is for the first day
                    if (recorded==0) and (hist_marker==1):
                        year_rec[year][dt.datetime(year,mon,day,hour)-dt.timedelta(hours=1)]=[re_cord_lon, re_cord_lat]
    
                    # This condition is for the remaining days
                    if (len(latss)>0) and (hist_marker==1) :
                        year_rec[year][dt.datetime(year,mon,day,hour)]=[lonss, latss]
                        recorded=1
                
                    elif ((len(lats)>0) and (hist_marker==0)) :
                        s10_hist=s10
                        u_hist=u
                        v_hist=v
                        lons_hist=lons
                        lats_hist=lats
                        recorded=0
                                                          
                    # Record lat lon values for the next day 
                    rec_cord=[]
                    for i in range(len(lats)):
                        rec_cord.append([lats[i],lons[i]]) '''
                        
        fields=['year','month','day','hour','longitude','latitude','maximum_wind','MW_lon','MW_lat','MW_radius']
    
        rows=[]
        print('creating csv file')
        #for year in tg_tqdm(year_rec, '1582179762:AAGbwgXCsNhQUWT35qGNGkTorAs-xa9Cwfk', '1250857525',desc='creating csv file') :
        for time in year_rec[year]:
            for i in range(len(year_rec[year][time][0])):
                line=[]
                line.append(time.year)
                line.append(time.month)
                line.append(time.day)
                line.append(time.hour)
                line.append(year_rec[year][time][0][i])
                line.append(year_rec[year][time][1][i])
                line=line+max_wind(year_rec[year][time][0][i],year_rec[year][time][1][i],time.year,time.month,time.day,time.hour,res)
                rows.append(line)
                                                                                                                
        import csv
        filename = '/home/picard/git_repo/results/mov_records_'+str(year)+'.csv'
        with open(filename, 'w') as csvfile: 
            csvwriter = csv.writer(csvfile)  
            csvwriter.writerow(fields)
            csvwriter.writerows(rows)
        csvfile.close()
