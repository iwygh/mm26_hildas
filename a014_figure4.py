#%%
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import time
#%%
t0 = time.time()
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
zinvar = np.sqrt(1-xinvar**2-yinvar**2)
idegplanets = np.loadtxt('b003_planets_ideg_'+yrsstr+'.csv',delimiter=',')
Wdegplanets = np.loadtxt('b003_planets_nodedeg_'+yrsstr+'.csv',delimiter=',')
nt = len(tyrsvec)
sini = np.sin(np.radians(idegplanets))
cosW = np.cos(np.radians(Wdegplanets))
sinW = np.sin(np.radians(Wdegplanets))
qplanets = sini * cosW
pplanets = sini * sinW
xJ = qplanets[0,:]
yJ = pplanets[0,:]
zJ = np.sqrt(1-xJ**2-yJ**2)
#%%
nseeds = 4
nclones = 500
seeds = ['Schubart','Ismene','Hilda','Potomac']
seednumbers = ['1911','190','153','1345']
short_labels = ['S','BPHSFG','H','P']
colors = ['blue','green','red','goldenrod']
input_conditions = [['Schubart'],['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],\
                    ['Hilda'],['Potomac']]
iobj_list = []
adjusts = [[-1.3,-1.3,-1.3],[-1.3,-1.3,-1.3],[-1.4,-1.4,-1.4],[+1.6,+1.6,+1.6]]
xcc_vmf_mat = np.zeros((nt,nseeds))
ycc_vmf_mat = np.zeros((nt,nseeds))
zcc_vmf_mat = np.zeros((nt,nseeds))
xL_mat = np.zeros((nt,nseeds))
yL_mat = np.zeros((nt,nseeds))
zL_mat = np.zeros((nt,nseeds))
x0_mat = np.zeros((nt,nseeds))
y0_mat = np.zeros((nt,nseeds))
z0_mat = np.zeros((nt,nseeds))
for iseed in range(nseeds):
    short_label = short_labels[iseed]
    seed = seeds[iseed]
    seednumber = seednumbers[iseed]
    df = pd.read_csv('b008_'+short_label+'_statistics_vmf_'+yrsstr+'.csv')
    xcc = np.array(df['xcc'].to_list())
    ycc = np.array(df['ycc'].to_list())
    zcc = np.sqrt(1-xcc**2-ycc**2)
    xcc_vmf_mat[:,iseed] = xcc
    ycc_vmf_mat[:,iseed] = ycc
    zcc_vmf_mat[:,iseed] = zcc
    dfl_icon = pd.read_csv('b007_laplaceplane_'+short_label+'_'+yrsstr+'.csv')
    dflx = np.array(dfl_icon['laplacex'].to_list())
    dfly = np.array(dfl_icon['laplacey'].to_list())
    dflz = np.sqrt(1-dflx**2-dfly**2)
    xL_mat[:,iseed] = dflx
    yL_mat[:,iseed] = dfly
    zL_mat[:,iseed] = dflz
    input_groups = input_conditions[iseed]
    ideg_here = np.loadtxt('b004_astorb_ideg_'+seed+'_'+yrsstr+'.csv',delimiter=',')
    Wdeg_here = np.loadtxt('b004_astorb_nodedeg_'+seed+'_'+yrsstr+'.csv',delimiter=',')
    irad_here = np.radians(ideg_here)
    Wrad_here = np.radians(Wdeg_here)
    q_here = np.sin(irad_here)*np.cos(Wrad_here)
    p_here = np.sin(irad_here)*np.sin(Wrad_here)
    s_here = np.cos(irad_here)
    x0_mat[:,iseed] = q_here
    y0_mat[:,iseed] = p_here
    z0_mat[:,iseed] = s_here
#%%
xinvar = xinvar[0]
yinvar = yinvar[0]
zinvar = zinvar[0]
nrefs = 3
nice_titles = ['Laplace Plane','Jupiter Plane','Invariable Plane']
short_titles = ['xL','xJ','xinvar']
xrefs = [xL_mat,xJ,xinvar]
yrefs = [yL_mat,yJ,yinvar]
zrefs = [zL_mat,zJ,zinvar]
svecs = np.zeros((nrefs,nseeds,nt))
fig2 = plt.figure(figsize=(12,4)) # (width,height)
nrows = 1
ncols = 3
ax2 = [fig2.add_subplot(nrows,ncols,i+1) for i in range(nrefs)]
plt.rcParams.update({'font.size':10})
for iax in range(nrefs):
    a = ax2[iax]
    a.tick_params(axis='x',direction='in')
    a.tick_params(axis='y',direction='in')
    a.set_xlim([0,2])
    a.set_ylim([0,15])
fig2.subplots_adjust(wspace=0.1, hspace=0.2)
tvec_yrs = np.linspace(0,2e6,10001,endpoint=True)
for iref in range(nrefs):
    nice_title = nice_titles[iref]
    short_title = short_titles[iref]
    xlim = [0,2]
    ylim = [0,15]
    xref = xrefs[iref]
    yref = yrefs[iref]
    zref = zrefs[iref]
    ax2[iref].set_title('Inclination to ' + nice_title)
    for iseed in range(nseeds):
        seednumber = seednumbers[iseed]
        adjust = adjusts[iseed][iref]
        seed = seeds[iseed]
        x0 = x0_mat[:,iseed]
        y0 = y0_mat[:,iseed]
        z0 = z0_mat[:,iseed]
        if iref in [0]:
            xref_here = xref[:,iseed]
            yref_here = yref[:,iseed]
            zref_here = zref[:,iseed]
        else: 
            xref_here = xref
            yref_here = yref
            zref_here = zref
        dotproduct = x0*xref_here + y0*yref_here + z0*zref_here
        svec = np.degrees(np.arccos(dotproduct))
        svecs[iref,iseed,:] = svec
        ax2[iref].plot(tyrsvec/1e6,svec,linewidth=0.2,c=colors[iseed])
        ds = np.nanmax(svec)-np.nanmin(svec)
        dm = np.nanmean(svec)
        ds2 = np.nanmax(np.abs(svec-np.nanmean(svec)))
        ax2[iref].axhline(np.nanmean(svec),linewidth=1,c='black')
        print(short_title, seed)
        print('mean i_free = ',str(np.round(dm,3)),' deg')
        print('Delta i_free = ',str(np.round(ds2,3)),' deg')
        print('')
        ax2[iref].text(1.9,np.nanmean(svec)+adjust,seednumber+' '+seed,horizontalalignment='right',\
                verticalalignment='center',bbox=dict(facecolor='white',alpha=0))
ax2[1].set_yticklabels([])
ax2[2].set_yticklabels([])
ax2[0].set_xlabel('Time (Myr)')
ax2[1].set_xlabel('Time (Myr)')
ax2[2].set_xlabel('Time (Myr)')
ax2[0].set_ylabel(r"$i_{rel}$"+' (degrees)')
fig2.savefig('figure4.png',dpi=400)
t1 = time.time()
print('runtime = ',str(np.round((t1-t0)/60,3)),' minutes')