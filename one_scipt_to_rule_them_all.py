import os
import pandas as pd

for year in range(1982,2021):
    imd_data=pd.read_csv('/home/picard/sas/mov/data/ecmwf/imd_files/'+str(year)+'.csv')
    if len(imd_data)>0:
        df=pd.read_csv('/home/picard/git_repo/results/mov_records_'+str(year)+'.csv')
        df=df.drop(columns=df.columns[df.columns.str.contains('error')])
        df.to_csv('/home/picard/git_repo/results/mov_records_'+str(year)+'.csv',index=False)

os.system('python3 /home/picard/git_repo/synthetic_vortex.py')
os.system('python3 /home/picard/git_repo/create_data_list.py')
os.system('python3 /home/picard/git_repo/plot_error.py')