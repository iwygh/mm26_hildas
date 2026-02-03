#%%
def vmf_fun(xvec,yvec,zvec,probmass):
    import numpy as np
    nobj = len(xvec)
    Sx = np.sum(xvec)
    Sy = np.sum(yvec)
    Sz = np.sum(zvec)
    R = np.sqrt(Sx**2+Sy**2+Sz**2)
    xcc = Sx/R
    ycc = Sy/R
    zcc = Sz/R
    Rbar = R/nobj
    dsum = 0
    for iobj in range(nobj):
        dsum = dsum + (xvec[iobj]*xcc+yvec[iobj]*ycc+zvec[iobj]*zcc)**2
    d = 1 - 1/nobj * dsum
    sigmahat = np.sqrt(d/(nobj*Rbar**2))
    A = 1 - probmass
    angledeg = np.degrees(np.arcsin(sigmahat*np.sqrt(-np.log(A))))
    return xcc,ycc,zcc,angledeg
#%%
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import time
#%%
jd = 2460796.5 # May 1, 2025 00:00:00
time_yrs  = int(2e6)
tstep_yrs = int(2e2)
dt_yrs = 0.2
tyrs_str = '2e6yr'
tstepyrs_str = '2e2yr'
dtyrs_str = '0.2yr'
yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
tyrsvec = np.loadtxt('b003_tyrsvec_' + yrsstr + '.csv',delimiter=',')
nt = len(tyrsvec)
dfa = pd.read_csv('b006_astorb_labels_Hmax16.3.csv')
n_astorb = dfa.shape[0]
Hmax = 16.3
ideg_mat = np.loadtxt('b004_astorb_ideg_0_'+yrsstr+'.csv',delimiter=',')
Wdeg_mat = np.loadtxt('b004_astorb_nodedeg_0_'+yrsstr+'.csv',delimiter=',')
output_labels = ['BPHSFG','H','S','P']   
plot_labels = ['Hilda group','Hilda family','Schubart family','Potomac family']
input_conditions = [['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],\
                    ['Hilda'],['Schubart'],['Potomac']]
ncons = len(input_conditions)
nreps = 100000
probmasses = [0.68,0.95,0.997]
t0 = time.time()
xcc_vecs = np.zeros((ncons,nreps))
ycc_vecs = np.zeros((ncons,nreps))
zcc_vecs = np.zeros((ncons,nreps))
angle_68_deg_vecs = np.zeros((ncons,nreps))
angle_95_deg_vecs = np.zeros((ncons,nreps))
angle_997_deg_vecs = np.zeros((ncons,nreps))
Delta_deg_vecs = np.zeros((ncons,nreps))
fig = plt.figure(figsize=(7,6)) # (width,height)
nrows = 2
ncols = 2
nplots = nrows * ncols
ax = [fig.add_subplot(nrows,ncols,i+1) for i in range(nrows*ncols)]
plt.rcParams.update({'font.size':10})
for iax in range(nplots):
    a = ax[iax]
    a.tick_params(axis='x',direction='in')
    a.tick_params(axis='y',direction='in')
subplottitles = plot_labels
for iax in range(nplots):
    ax[iax].text(0.9,0.9,subplottitles[iax],horizontalalignment='right',verticalalignment='top',\
               transform=ax[iax].transAxes,bbox=dict(facecolor='white',alpha=1))
    ax[iax].set_xlabel('Relative inclination (deg)')
    ax[iax].set_ylabel('Count')
fig.subplots_adjust(wspace=0.3, hspace=0.3)
for icon in range(ncons):
    output_label = output_labels[icon]
    plot_label = plot_labels[icon]
    input_groups = input_conditions[icon]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) :
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(output_label,nobj)
    irad_mat = np.radians(ideg_mat[astorb_indices])
    Wrad_mat = np.radians(Wdeg_mat[astorb_indices])
    xvec = np.sin(irad_mat)*np.cos(Wrad_mat)
    yvec = np.sin(irad_mat)*np.sin(Wrad_mat)
    zvec = np.cos(irad_mat)
    xcc,ycc,zcc,angle1deg = vmf_fun(xvec,yvec,zvec,probmasses[0])
    xrel = xvec - xcc
    yrel = yvec - ycc
    irel = np.arcsin(np.sqrt(xrel**2 + yrel**2))
    ireldeg = np.degrees(irel)
    ax[icon].hist(ireldeg,bins=20,edgecolor='black',color='skyblue',alpha=0.8)
    angle_68_deg_vec = []
    angle_95_deg_vec = []
    angle_997_deg_vec = []
    Delta_deg_vec = []
    xcc_vec = []
    ycc_vec = []
    zcc_vec = []
    for irep in range(nreps):
        irel_irep = np.random.choice(irel,size=nobj,replace=True)
        Wrel_irep = np.random.uniform(low=0,high=2*np.pi,size=nobj)
        xvec_irep = np.sin(irel_irep)*np.cos(Wrel_irep)
        yvec_irep = np.sin(irel_irep)*np.sin(Wrel_irep)
        zvec_irep = np.cos(irel_irep)
        xcc_irep,ycc_irep,zcc_irep,angle1deg_irep = vmf_fun(xvec_irep,yvec_irep,zvec_irep,probmasses[0])
        xcc_irep,ycc_irep,zcc_irep,angle2deg_irep = vmf_fun(xvec_irep,yvec_irep,zvec_irep,probmasses[1])
        xcc_irep,ycc_irep,zcc_irep,angle3deg_irep = vmf_fun(xvec_irep,yvec_irep,zvec_irep,probmasses[2])
        xcc_vec.append(xcc_irep)
        ycc_vec.append(ycc_irep)
        zcc_vec.append(zcc_irep)
        angle_68_deg_vec.append(angle1deg_irep)
        angle_95_deg_vec.append(angle2deg_irep)
        angle_997_deg_vec.append(angle3deg_irep)
        Delta_deg_vec.append(np.degrees(np.arccos(zcc_irep)))
        t1 = time.time()
    angle_68_deg_vec = np.array(angle_68_deg_vec)
    angle_95_deg_vec = np.array(angle_95_deg_vec)
    angle_997_deg_vec = np.array(angle_997_deg_vec)
    Delta_deg_vec = np.array(Delta_deg_vec)
    xcc_vec = np.array(xcc_vec)
    ycc_vec = np.array(ycc_vec)
    zcc_vec = np.array(zcc_vec)
    xcc_vecs[icon,:] = xcc_vec
    ycc_vecs[icon,:] = ycc_vec
    zcc_vecs[icon,:] = zcc_vec
    angle_68_deg_vecs[icon,:] = angle_68_deg_vec
    angle_95_deg_vecs[icon,:] = angle_95_deg_vec
    angle_997_deg_vecs[icon,:] = angle_997_deg_vec
    Delta_deg_vecs[icon,:] = Delta_deg_vec
plt.savefig('figure1.png',dpi=400)
print(np.sum(Delta_deg_vecs[0,:]<angle_68_deg_vecs[0,:])/nreps)
print(np.sum(Delta_deg_vecs[1,:]<angle_68_deg_vecs[1,:])/nreps)
print(np.sum(Delta_deg_vecs[2,:]<angle_68_deg_vecs[2,:])/nreps)
print(np.sum(Delta_deg_vecs[3,:]<angle_68_deg_vecs[3,:])/nreps)
print('')
print(np.sum(Delta_deg_vecs[0,:]<angle_95_deg_vecs[0,:])/nreps)
print(np.sum(Delta_deg_vecs[1,:]<angle_95_deg_vecs[1,:])/nreps)
print(np.sum(Delta_deg_vecs[2,:]<angle_95_deg_vecs[2,:])/nreps)
print(np.sum(Delta_deg_vecs[3,:]<angle_95_deg_vecs[3,:])/nreps)
print('')
print(np.sum(Delta_deg_vecs[0,:]<angle_997_deg_vecs[0,:])/nreps)
print(np.sum(Delta_deg_vecs[1,:]<angle_997_deg_vecs[1,:])/nreps)
print(np.sum(Delta_deg_vecs[2,:]<angle_997_deg_vecs[2,:])/nreps)
print(np.sum(Delta_deg_vecs[3,:]<angle_997_deg_vecs[3,:])/nreps)
    