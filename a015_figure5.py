#%%
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
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
xinvar = np.loadtxt('b003_planets_xinvar_'+yrsstr+'.csv',delimiter=',')
yinvar = np.loadtxt('b003_planets_yinvar_'+yrsstr+'.csv',delimiter=',')
idegplanets = np.loadtxt('b003_planets_ideg_'+yrsstr+'.csv',delimiter=',')
Wdegplanets = np.loadtxt('b003_planets_nodedeg_'+yrsstr+'.csv',delimiter=',')
sini = np.sin(np.radians(idegplanets))
cosW = np.cos(np.radians(Wdegplanets))
sinW = np.sin(np.radians(Wdegplanets))
qplanets = sini * cosW
pplanets = sini * sinW
xJ = qplanets[0,:]
yJ = pplanets[0,:]
output_labels = ['P','H','BPHSFG','S']
input_conditions = [['Potomac'],['Hilda'],['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['Schubart']]
clone_labels = ['Potomac','Hilda','Ismene','Schubart']
fig = plt.figure(figsize=(5,8)) # (width,height)
nrows = 4
ncols = 2
ax = [fig.add_subplot(nrows,ncols,i+1) for i in range(nrows*ncols)]
plt.rcParams.update({'font.size':6})
for iax in [0,1,2,3,4,5,6,7]:
    a = ax[iax]
    a.tick_params(axis='x',direction='in')
    a.tick_params(axis='y',direction='in')
fig.subplots_adjust(wspace=0.05, hspace=0.05)
subplottitles = ['Potomac family',   'Potomac clones',\
                 'Hilda family',     'Hilda clones',\
                 'Hilda group',  'Ismene clones',\
                 'Schubart family','Schubart clones']
for iax in [0,1,2,3,4,5,6,7]:
    ax[iax].text(0.06,0.9,subplottitles[iax],horizontalalignment='left',verticalalignment='top',\
               transform=ax[iax].transAxes,bbox=dict(facecolor='white',alpha=1))
    ax[iax].grid('on')
for iax in [0,2,4,6]:
    averaging_window = 250
    output_label = output_labels[int(iax/2)]
    input_condition = input_conditions[int(iax/2)]
    dfvmf = pd.read_csv('b008_'+output_label+'_statistics_vmf_'+yrsstr+'.csv')
    xcc = np.array(dfvmf['xcc'].to_list())
    ycc = np.array(dfvmf['ycc'].to_list())
    angle95deg = np.array(dfvmf['angle_95_deg'].to_list())
    dfl_icon = pd.read_csv('b007_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
    dflx = dfl_icon['laplacex'].to_list()
    dfly = dfl_icon['laplacey'].to_list()
    deltadeg_laplace_here = np.degrees(np.arcsin(np.sqrt((xcc-dflx)**2+(ycc-dfly)**2)))
    deltadeg_jupiter_here = np.degrees(np.arcsin(np.sqrt((xcc-xJ)**2+(ycc-yJ)**2)))
    deltadeg_invariable_here = np.degrees(np.arcsin(np.sqrt((xcc-xinvar)**2+(ycc-yinvar)**2)))
    deltadeg_laplace_here_av = []
    deltadeg_jupiter_here_av = []
    deltadeg_invariable_here_av = []
    angle95deg_av = []
    for it in range(nt-averaging_window):
        ddlh = np.mean(deltadeg_laplace_here[it:it+averaging_window])
        ddjh = np.mean(deltadeg_jupiter_here[it:it+averaging_window])
        ddih = np.mean(deltadeg_invariable_here[it:it+averaging_window])
        dd9h = np.mean(angle95deg[it:it+averaging_window])
        deltadeg_laplace_here_av.append(ddlh)
        deltadeg_jupiter_here_av.append(ddjh)
        deltadeg_invariable_here_av.append(ddih)
        angle95deg_av.append(dd9h)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_laplace_here_av,c='orange',lw=0.25)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_jupiter_here_av,c='black',lw=0.25)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_invariable_here_av,c='purple',lw=0.25)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,angle95deg_av,c='red',lw=0.25)
    ax[iax].grid('on')
    ax[iax].set_xlim([0,2])
    ax[iax].set_ylabel(r"$\Delta i$"+' (degrees)')
for iax in [1,3,5,7]:
    averaging_window = 250
    clone_label = clone_labels[int((iax-1)/2)]
    output_label = output_labels[int((iax-1)/2)]
    input_condition = input_conditions[int((iax-1)/2)]
    dfcf = pd.read_csv('b009_'+clone_label+'_clones_statistics_cf_Hmax16.3_'+yrsstr+'.csv')
    xcc = np.array(dfcf['xcc'].to_list())
    ycc = np.array(dfcf['ycc'].to_list())
    angle95deg = np.array(dfcf['sigs_95_deg'].to_list())
    dfl_icon = pd.read_csv('b007_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
    dflx = dfl_icon['laplacex'].to_list()
    dfly = dfl_icon['laplacey'].to_list()
    deltadeg_laplace_here = np.degrees(np.arcsin(np.sqrt((xcc-dflx)**2+(ycc-dfly)**2)))
    deltadeg_jupiter_here = np.degrees(np.arcsin(np.sqrt((xcc-xJ)**2+(ycc-yJ)**2)))
    deltadeg_invariable_here = np.degrees(np.arcsin(np.sqrt((xcc-xinvar)**2+(ycc-yinvar)**2)))
    deltadeg_laplace_here_av = []
    deltadeg_jupiter_here_av = []
    deltadeg_invariable_here_av = []
    angle95deg_av = []
    for it in range(nt-averaging_window):
        ddlh = np.mean(deltadeg_laplace_here[it:it+averaging_window])
        ddjh = np.mean(deltadeg_jupiter_here[it:it+averaging_window])
        ddih = np.mean(deltadeg_invariable_here[it:it+averaging_window])
        dd9h = np.mean(angle95deg[it:it+averaging_window])
        deltadeg_laplace_here_av.append(ddlh)
        deltadeg_jupiter_here_av.append(ddjh)
        deltadeg_invariable_here_av.append(ddih)
        angle95deg_av.append(dd9h)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_laplace_here_av,c='orange',lw=0.25)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_jupiter_here_av,c='black',lw=0.25)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,deltadeg_invariable_here_av,c='purple',lw=0.25)
    ax[iax].plot(tyrsvec[:nt-averaging_window]/1e6,angle95deg_av,c='red',lw=0.25)
    ax[iax].grid('on')
ax[0].set_ylim([0,1])
ax[2].set_ylim([0,1])
ax[4].set_ylim([0,1])
ax[6].set_ylim([0,1])
ax[1].set_ylim([0,1])
ax[3].set_ylim([0,1])
ax[5].set_ylim([0,1])
ax[7].set_ylim([0,1])
ax[0].set_xlim([0,2])
ax[2].set_xlim([0,2])
ax[4].set_xlim([0,2])
ax[6].set_xlim([0,2])
ax[1].set_xlim([0,2])
ax[3].set_xlim([0,2])
ax[5].set_xlim([0,2])
ax[7].set_xlim([0,2])
ax[6].set_xticks([0,0.5,1,1.5])
ax[7].set_xticks([0,0.5,1,1.5])
ax[1].set_yticks([0,0.25,0.5,0.75])
ax[3].set_yticks([0,0.25,0.5,0.75])
ax[5].set_yticks([0,0.25,0.5,0.75])
ax[7].set_yticks([0,0.25,0.5,0.75])
ax[0].set_yticks([0,0.25,0.5,0.75])
ax[2].set_yticks([0,0.25,0.5,0.75])
ax[4].set_yticks([0,0.25,0.5,0.75])
ax[6].set_yticks([0,0.25,0.5,0.75])
ax[1].set_yticklabels([])
ax[3].set_yticklabels([])
ax[5].set_yticklabels([])
ax[7].set_yticklabels([])
for iax in range(6):
    ax[iax].set_xticklabels([])
ax[6].set_xlabel('Time, Myr')
ax[7].set_xlabel('Time, Myr')
outfile_plot = 'figure5.png'
plt.savefig(outfile_plot,dpi=400)