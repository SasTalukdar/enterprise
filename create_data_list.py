import pandas as pd

df_list=[]
for year in range(1982,2021):
    imd_data=pd.read_csv('/home/picard/sas/mov/data/ecmwf/imd_files/'+str(year)+'.csv')
    if len(imd_data)>0:
        df=pd.read_csv('/home/picard/git_repo/results/mov_records_'+str(year)+'.csv')
        df_list.append(df)

df=df_list[0]
for i in range(1,len(df_list)):
    df=pd.concat([df,df_list[i]],axis=0)

df=df.reset_index()
df=df.drop(columns='index')

df.to_csv('/home/picard/git_repo/results/mov_data.csv',index=False)