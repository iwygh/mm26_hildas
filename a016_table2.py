#%%
import numpy as np
import pandas as pd
jd = 2460796.5 # May 1, 2025 00:00:00
time_yrs  = int(2e6)
tstep_yrs = int(2e2)
dt_yrs = 0.2
tyrs_str = '2e6yr'
tstepyrs_str = '2e2yr'
dtyrs_str = '0.2yr'
yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
#%%
output_labels = ['BPHSFG','H','S','P']
labels_files = ['b006_all_Hmax16.3.csv','b006_hilda_Hmax16.3.csv','b006_schubart_Hmax16.3.csv','b006_potomac_Hmax16.3.csv']
nlabels = len(output_labels)
for ilabel in range(nlabels):
    labels_file = labels_files[ilabel]
    label = output_labels[ilabel]
    df = pd.read_csv(labels_file)
    nobj = df.shape[0]
    vmf_file = 'b008_'+label+'_statistics_vmf_'+yrsstr+'.csv'
    df = pd.read_csv(vmf_file)
    xcc = df['xcc'][0]
    ycc = df['ycc'][0]
    zcc = np.sqrt(1-xcc**2-ycc**2)
    i0rad = np.arccos(zcc)
    i0deg = np.degrees(i0rad)
    sini = np.sin(i0rad)
    W0rad = np.arctan2(ycc/sini,xcc/sini)
    W0deg = np.degrees(W0rad)
    phi95deg = df['angle_95_deg'][0]
    print(label,nobj,np.round(i0deg,2),np.round(W0deg,2),np.round(phi95deg,2))
