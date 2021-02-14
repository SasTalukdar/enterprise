import numpy as np

def normalize(u,v):
    mod=np.sqrt(u*u+v*v)
    u=u/mod
    v=v/mod
    return [u,v]

def t(m,d,h):
    if m==5:
        t1=0
    elif m==6:
        t1=31
    t=t1*24+(d-1)*24+h+1
    return t

def dist(i,j):
    d=((i-mid_i)**2+(j-mid_j)**2)**0.5
    d=d*res
    return d

def error(u,v,u_ref,v_ref):
    [u,v]=normalize(u,v)
    e=0
    for i in range(np.shape(u)[0]):
        for j in range(np.shape(u)[1]):
            if dist(i,j)<=radius:
                e=e+((u[i][j]-u_ref[i][j])**2+(v[i][j]-v_ref[i][j])**2)**0.5
    e=e/((i+1)*(j+1))
    return e

for radius in range(100,800,50):
    mid=int((radius/111.1)*4)
    u_ref=np.zeros((2*mid+1,2*mid+1))
    v_ref=np.zeros((2*mid+1,2*mid+1))
    mid_i=mid
    mid_j=mid
    #r_m=(389.489/111.1)*4
    r_m=1
    #v_ref_m=9.89
    v_ref_m=1
    b=1
    for i in range(len(u_ref)):
        for j in range(len(u_ref)):
            i_o=i-mid_i
            j_o=j-mid_j
            r=np.sqrt(i_o**2+j_o**2)
            v_ref_t=v_ref_m*(r/r_m)*np.exp((1/b)*(1-(r/r_m)**b))
            u_ref[i][j]=v_ref_t*i_o/r
            v_ref[i][j]=v_ref_t*j_o/r
    
    [u_ref,v_ref]=normalize(u_ref, v_ref)
    u_ref=np.nan_to_num(u_ref)
    v_ref=np.nan_to_num(v_ref)
    u_ref=np.flip(u_ref)
    #v_ref=np.flip(v_ref)
    
    #x,y=np.meshgrid(range(29),range(29))
    #y=np.flip(y)
    
    #u_ref=u_ref/mod
    #v_ref=v_ref/mod
    #s=np.sqrt(u_ref*u_ref+v_ref*v_ref)
    
    #plt.streamplot(x, y, u_ref, v_ref)
    #plt.qu_refiv_refer(x,y,u_ref,v_ref,headwidth=5,headlength=10)
    
    import pandas as pd
    from py3grads import Grads
    ga=Grads(verbose=False)
    
    res=27.75
    run=6
    for year in range(1982,2021):
        imd_data=pd.read_csv('/home/picard/sas/mov/data/ecmwf/imd_files/'+str(year)+'.csv')
        if len(imd_data)>0:
            df=pd.read_csv('/home/picard/git_repo/results/mov_records_'+str(year)+'.csv')
            err=[]
            for i in df.index:
                t2=t(df.month[i],df.day[i],df.hour[i])
                u=[]
                v=[]
                ga('reinit')
                ga('sdfopen /home/picard/sas/mov/data/ecmwf/u_v_gp/era5_hourly_u_v_gp_850_'+str(year)+'.nc')
                ga('set lat '+str(df.latitude[i]-mid/4)+' '+str(df.latitude[i]+mid/4))
                ga('set lon '+str(df.longitude[i]-mid/4)+' '+str(df.longitude[i]+mid/4))
                for time in range(t2-run,t2+run+1):
                    ga('set t '+str(time))
                    u.append(ga.exp('u'))
                    v.append(ga.exp('v'))
                u=np.array(u)
                u=np.nanmean(u,axis=0)
                v=np.array(v)
                v=np.nanmean(v,axis=0)
                err.append(error(u,v,u_ref,v_ref))
            df['error_index_'+str(radius)]=err
            
            df.to_csv('/home/picard/git_repo/results/mov_records_'+str(year)+'.csv',index=False)