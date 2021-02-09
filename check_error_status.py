import pandas as pd

for year in range(1982,2021):
    imd_data=pd.read_csv('/home/picard/sas/mov/data/ecmwf/imd_files/'+str(year)+'.csv')
    if len(imd_data)>0:
        df_ref=pd.read_csv('/home/picard/sas/mov/data/ecmwf/mean_results/run_mean_6_hour_mov_records_'+str(year)+'.csv')
        ref_list=[]
        for i,row in df_ref.iterrows():
            ref_list.append([row.month,row.day,row.hour,row.latitude,row.longitude])
        
        df=pd.read_csv('/home/picard/git_repo/results/mov_records_'+str(year)+'.csv')
        df['check']=0
        for i,row in df.iterrows():
            if [row.month,row.day,row.hour,row.latitude,row.longitude] in ref_list:
                df['check'][i]=1
        df.to_csv('/home/picard/git_repo/results/mov_records_'+str(year)+'.csv',index=False)