import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df=pd.read_csv('results/mov_data.csv')

for radius in range(100,800,50):
    th_value=[]
    pod=[]
    far=[]
    tp=len(df[df.check==1])
    tn=len(df[df.check==0])
    for th in np.arange(np.nanmin(df['error_index_'+str(radius)])+0.01,1,0.01):
        th_value.append(th)
        p=0
        fn=0
        all_p=0
        for i,row in df.iterrows():
            if row['error_index_'+str(radius)]<=th and row.check==1:
                p=p+1
            elif row['error_index_'+str(radius)]<=th and row.check==0:
                fn=fn+1
        pod.append(p/tp)
        far.append(fn/tn)
    plt.plot(th_value,pod,label='POD')
    plt.plot(th_value,far,label='FAR')
    plt.legend()
    plt.xlabel('Threshold Value for Error Index')
    plt.xlim(0,1)
    plt.xticks(np.arange(0,1.05,0.05),rotation=45)
    plt.yticks(np.arange(0,1.05,0.05))
    plt.grid()
    plt.title('Analysis based on Error Index at '+str(radius)+' km')
    plt.tight_layout()
    plt.savefig('/home/picard/git_repo/error_index_analysis/error_index '+str(radius)+'.png',dpi=300)
    plt.clf()
