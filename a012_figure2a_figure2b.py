#%%
def plot_lines_and_spokes(ax):
    import numpy as np
    th = np.linspace(start=0,stop=2*np.pi,num=100,endpoint=True)
    costh = np.cos(th)
    sinth=  np.sin(th)
    for angle in [90,80,70,60,50,40,30,20,10,0]: # 10-deg latitude lines
        ax.plot(np.cos(np.radians(angle))*costh,np.cos(np.radians(angle))*sinth,color='gray',\
                linestyle='-',linewidth=0.25)
    for integer in [0,1,2,3,4,5,6,7,8,9,10,11]: # 30-deg longitude lines
        ax.plot([0,np.cos(np.radians(integer*30))],[0,np.sin(np.radians(integer*30))],color='gray',\
                linestyle='-',linewidth=0.25)
    # x and y axis lines
    ax.axhline(color='black',linestyle='-',linewidth=1)
    ax.axvline(color='black',linestyle='-',linewidth=1)
    return ax
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
aau_mat_planets = np.loadtxt('b003_planets_aau_'+yrsstr+'.csv',delimiter=',')
tyrsvec = np.loadtxt('b003_tyrsvec_' + yrsstr + '.csv',delimiter=',')
nt = len(tyrsvec)
dfa = pd.read_csv('b006_astorb_labels_Hmax16.3.csv')
n_astorb = dfa.shape[0]
probmasses = [0.68,0.95,0.997]
Hmax = 16.3
aau_astorb = np.loadtxt('b004_astorb_aau_0_'+yrsstr+'.csv',delimiter=',')
e_astorb = np.loadtxt('b004_astorb_e_0_'+yrsstr+'.csv',delimiter=',')
ideg_astorb = np.loadtxt('b004_astorb_ideg_0_'+yrsstr+'.csv',delimiter=',')
Wdeg_astorb = np.loadtxt('b004_astorb_nodedeg_0_'+yrsstr+'.csv',delimiter=',')
sini = np.sin(np.radians(ideg_astorb))
cosW = np.cos(np.radians(Wdeg_astorb))
sinW = np.sin(np.radians(Wdeg_astorb))
q_astorb = sini*cosW
p_astorb = sini*sinW
output_labels = ['P','H','S','B','BFG','BPHS','BPHSFG']
input_conditions = [['Potomac'],['Hilda'],['Schubart'],['Background'],\
                    ['Background','Francette','Guinevere'],['Background','Potomac','Hilda','Schubart'],\
                    ['Background','Potomac','Hilda','Schubart','Francette','Guinevere']]
colors = ['goldenrod','red','blue','green','green','black','black']   
distinct_groups = ['P','H','S','B','BFG','BPHS','BPHSFG']
conditions_for_distinct_groups = [[['Potomac'],['librating']],\
    [['Hilda'],['librating']],\
    [['Schubart'],['librating']],\
    [['Background'],['librating']],[['Background'],['maybe']],[['Background'],['circulating']],\
    [['Background','Francette','Guinevere'],['librating']],[['Background','Francette','Guinevere'],['maybe']],[['Background','Francette','Guinevere'],['circulating']],\
    [['Background','Potomac','Hilda','Schubart'],['librating']],[['Background','Potomac','Hilda','Schubart'],['maybe']],[['Background','Potomac','Hilda','Schubart'],['circulating']],\
    [['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['librating']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['maybe']],[['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],['circulating']],\
    ] 
groups_in_conditions = [['P_L'],     ['H_L'],     ['S_L'],\
     ['B_L'],     ['B_L','B_M'],          ['B_L','B_M','B_C'],\
     ['BFG_L'],   ['BFG_L','BFG_M'],      ['BFG_L','BFG_M','BFG_C'],\
     ['BPHS_L'],  ['BPHS_L','BPHS_M'],    ['BPHS_L','BPHS_M','BPHS_C'],\
     ['BPHSFG_L'],['BPHSFG_L','BPHSFG_M'],['BPHSFG_L','BPHSFG_M','BPHSFG_C']] 
sizes = [5,20,40]
alphas = [0.1,0.3,1]
markers = ['.','s','o']
ncons = len(input_conditions)
#%%
outfile_plot = 'figure2b.png'
it = 0
labels_to_use = ['P','H','S','BFG']
nuse = len(labels_to_use)
xlim = [-0.25,0.25]
ylim = [-0.25,0.25]
fig = plt.figure(figsize=(7,7))
plt.rcParams['font.size'] = 18
ax = fig.add_subplot(111)
ax = plot_lines_and_spokes(ax)
for iuse in range(nuse):
    label_use = labels_to_use[iuse]
    indx = output_labels.index(label_use)
    input_groups = input_conditions[indx][0]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) :
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(label_use,nobj)
    aau_icon = aau_astorb[astorb_indices]
    e_icon = e_astorb[astorb_indices]
    ideg_icon = ideg_astorb[astorb_indices]
    q_icon = q_astorb[astorb_indices]
    p_icon = p_astorb[astorb_indices]
    x = q_icon
    y = p_icon
    ax.scatter(x,y,color=colors[indx],s=20,alpha=0.5,edgecolor='none',linewidth=0,marker='.')
plt.axis('equal')
ax.set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
ax.set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_box_aspect(1)
plt.tight_layout()
plt.savefig(outfile_plot,dpi=300)
plt.show()    
#%%
outfile_plot = 'figure2a.png'
it = 0
labels_to_use = ['BFG','P','H','S']
nuse = len(labels_to_use)
fig = plt.figure(figsize=(6,6))
plt.rcParams['font.size'] = 18
ax = fig.add_subplot(223)
for iuse in range(nuse):
    label_use = labels_to_use[iuse]
    indx = output_labels.index(label_use)
    input_groups = input_conditions[indx][0]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) :
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(label_use,nobj)
    aau_icon = aau_astorb[astorb_indices]
    e_icon = e_astorb[astorb_indices]
    ideg_icon = ideg_astorb[astorb_indices]
    x = aau_icon
    y = e_icon
    ax.scatter(x,y,color=colors[indx],s=10,alpha=0.5,edgecolor='none',linewidth=0,marker='.')
ax.set_xlabel('a')
ax.set_ylabel('e')
ax.tick_params(axis='x',direction='in')
ax.tick_params(axis='y',direction='in')
ax = fig.add_subplot(221)
for iuse in range(nuse):
    label_use = labels_to_use[iuse]
    indx = output_labels.index(label_use)
    input_groups = input_conditions[indx][0]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) :
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(label_use,nobj)
    aau_icon = aau_astorb[astorb_indices]
    e_icon = e_astorb[astorb_indices]
    ideg_icon = ideg_astorb[astorb_indices]
    x = aau_icon
    y = ideg_icon
    ax.scatter(x,y,color=colors[indx],s=10,alpha=0.5,edgecolor='none',linewidth=0,marker='.')
ax.set_xticklabels([])
ax.set_ylabel('i (deg)')
ax.tick_params(axis='x',direction='in')
ax.tick_params(axis='y',direction='in')
ax = fig.add_subplot(224)
for iuse in range(nuse):
    label_use = labels_to_use[iuse]
    indx = output_labels.index(label_use)
    input_groups = input_conditions[indx][0]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) :
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(label_use,nobj)
    aau_icon = aau_astorb[astorb_indices]
    e_icon = e_astorb[astorb_indices]
    ideg_icon = ideg_astorb[astorb_indices]
    x = e_icon
    y = ideg_icon
    ax.scatter(y,x,color=colors[indx],s=10,alpha=0.5,edgecolor='none',linewidth=0,marker='.')
ax.set_yticklabels([])
ax.tick_params(axis='x',direction='in')
ax.tick_params(axis='y',direction='in')
ax.set_xlabel('i (deg)')
fig.subplots_adjust(wspace=0, hspace=0)
plt.savefig(outfile_plot,dpi=300)
plt.show() 