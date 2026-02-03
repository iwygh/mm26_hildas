import numpy as np
import pandas as pd
#%%
jd = 2460796.5 # May 1, 2025 00:00:00
time_yrs  = int(2e6)
tstep_yrs = int(2e2)
dt_yrs = 0.2
tyrs_str = '2e6yr'
tstepyrs_str = '2e2yr'
dtyrs_str = '0.2yr'
yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
aau_astorb = np.loadtxt('b004_astorb_aau_'+yrsstr+'.csv',delimiter=',')
e_astorb = np.loadtxt('b004_astorb_e_'+yrsstr+'.csv',delimiter=',')
ideg_astorb = np.loadtxt('b004_astorb_ideg_'+yrsstr+'.csv',delimiter=',')
Wdeg_astorb = np.loadtxt('b004_astorb_nodedeg_'+yrsstr+'.csv',delimiter=',')
seeds = ['Schubart','Ismene','Hilda','Potomac']
nseeds = 4
dfa = pd.read_csv('b006_astorb_labels_Hmax16.3.csv')
idlist = dfa['idhere'].to_list()
n_astorb = dfa.shape[0]
nt = int(time_yrs/tstep_yrs)+1
#%%
for iseed in range(nseeds):
    seed = seeds[iseed]
    index = idlist.index(seed)
    print(seed,index)
    ideg_here = ideg_astorb[index,:]
    Wdeg_here = Wdeg_astorb[index,:]
    np.savetxt('b004_astorb_ideg_'+seed+'_'+yrsstr+'.csv',ideg_here,delimiter=',')
    np.savetxt('b004_astorb_nodedeg_'+seed+'_'+yrsstr+'.csv',Wdeg_here,delimiter=',')
it = 0
aau_astorb_0 = aau_astorb[:,it]
e_astorb_0 = e_astorb[:,it]
ideg_astorb_0 = ideg_astorb[:,it]
Wdeg_astorb_0 = Wdeg_astorb[:,it]
np.savetxt('b004_astorb_aau_0_'+yrsstr+'.csv',aau_astorb_0,delimiter=',')
np.savetxt('b004_astorb_e_0_'+yrsstr+'.csv',e_astorb_0,delimiter=',')
np.savetxt('b004_astorb_ideg_0_'+yrsstr+'.csv',ideg_astorb_0,delimiter=',')
np.savetxt('b004_astorb_nodedeg_0_'+yrsstr+'.csv',Wdeg_astorb_0,delimiter=',')

