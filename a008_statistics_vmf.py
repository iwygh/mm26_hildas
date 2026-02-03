#%%
def vmf_fun(xvec,yvec,zvec,probmass):
    import numpy as np
    nobj = len(xvec)
    assert len(xvec) == len(yvec) == len(zvec), 'x y z vecs not same length'
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
import astropy.stats.circstats as acstats
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
dfa = pd.read_csv('b006_astorb_labels_Hmax16.3.csv')
n_astorb = dfa.shape[0]
probmasses = [0.68,0.95,0.997]
Hmax = 16.3
ideg_mat = np.loadtxt('b004_astorb_ideg_'+yrsstr+'.csv',delimiter=',')
Wdeg_mat = np.loadtxt('b004_astorb_nodedeg_'+yrsstr+'.csv',delimiter=',')
output_labels = ['P','H','S','BPHSFG']   
input_conditions = [['Potomac'],['Hilda'],['Schubart'],\
                    ['Background','Potomac','Hilda','Schubart','Francette','Guinevere']]
ncons = len(input_conditions)
aau_mat_planets = np.loadtxt('b003_planets_aau_'+yrsstr+'.csv',delimiter=',')
xinvar_vec = np.loadtxt('b003_planets_xinvar_'+yrsstr+'.csv',delimiter=',')
yinvar_vec = np.loadtxt('b003_planets_yinvar_'+yrsstr+'.csv',delimiter=',')
idegplanets = np.loadtxt('b003_planets_ideg_'+yrsstr+'.csv',delimiter=',')
Wdegplanets = np.loadtxt('b003_planets_nodedeg_'+yrsstr+'.csv',delimiter=',')
sini = np.sin(np.radians(idegplanets))
cosW = np.cos(np.radians(Wdegplanets))
sinW = np.sin(np.radians(Wdegplanets))
qplanets = sini * cosW
pplanets = sini * sinW
xJ_vec = qplanets[0,:]
yJ_vec = pplanets[0,:]
for icon in range(ncons):
    output_label = output_labels[icon]
    input_groups = input_conditions[icon]
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) :
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    print(output_label,nobj)
    dfl = pd.read_csv('b007_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
    dflx_vec = dfl['laplacex'].to_list()
    dfly_vec = dfl['laplacey'].to_list()
    xcc_vec = []
    ycc_vec = []
    angle1deg_vec = []
    angle2deg_vec = []
    angle3deg_vec = []
    p_laplace_vec = []
    p_xcc_vec = []
    p_xJ_vec = []
    p_xI_vec = []
    irad_mat = np.radians(ideg_mat[astorb_indices,:])
    Wrad_mat = np.radians(Wdeg_mat[astorb_indices,:])
    xmat = np.sin(irad_mat)*np.cos(Wrad_mat)
    ymat = np.sin(irad_mat)*np.sin(Wrad_mat)
    zmat = np.cos(irad_mat)
    for it in range(nt):
        xvec = xmat[:,it]
        yvec = ymat[:,it]
        zvec = zmat[:,it]
        xcc,ycc,zcc,angle1deg = vmf_fun(xvec,yvec,zvec,probmasses[0])
        xcc,ycc,zcc,angle2deg = vmf_fun(xvec,yvec,zvec,probmasses[1])
        xcc,ycc,zcc,angle3deg = vmf_fun(xvec,yvec,zvec,probmasses[2])
        xcc_vec.append(xcc)
        ycc_vec.append(ycc)
        angle1deg_vec.append(angle1deg)
        angle2deg_vec.append(angle2deg)
        angle3deg_vec.append(angle3deg)
        dx = xvec - dflx_vec[it]
        dy = yvec - dfly_vec[it]
        dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
        dsini = np.sin(np.radians(dideg))
        dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
        pval = acstats.rayleightest(np.radians(dWdeg))
        p_laplace_vec.append(pval)
        dx = xvec - xcc
        dy = yvec - ycc
        dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
        dsini = np.sin(np.radians(dideg))
        dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
        pval = acstats.rayleightest(np.radians(dWdeg))
        p_xcc_vec.append(pval)
        dx = xvec - xJ_vec[it]
        dy = yvec - yJ_vec[it]
        dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
        dsini = np.sin(np.radians(dideg))
        dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
        pval = acstats.rayleightest(np.radians(dWdeg))
        p_xJ_vec.append(pval)
        dx = xvec - xinvar_vec[it]
        dy = yvec - yinvar_vec[it]
        dideg = np.degrees(np.arcsin(np.sqrt(dx**2+dy**2)))
        dsini = np.sin(np.radians(dideg))
        dWdeg = np.degrees(np.arctan2(dy/dsini,dx/dsini))
        pval = acstats.rayleightest(np.radians(dWdeg))
        p_xI_vec.append(pval)
    dictionary = {'xcc':xcc_vec,\
                  'ycc':ycc_vec,\
                  'angle_68_deg':angle1deg_vec,\
                  'angle_95_deg':angle2deg_vec,\
                  'angle_997_deg':angle3deg_vec,\
                  'pW_laplace_rayleightest':p_laplace_vec,\
                  'pW_xJ_rayleightest':p_xJ_vec,\
                  'pW_xinvar_rayleightest':p_xI_vec,\
                  'pW_xcc_rayleightest':p_xcc_vec,\
                      }
    dfd = pd.DataFrame.from_dict(dictionary)
    dfd.to_csv('b008_'+output_label+'_statistics_vmf_'+yrsstr+'.csv',index=None)
    plt.plot(tyrsvec/1e6,p_laplace_vec)
    plt.axhline(y=0.05,color='red')
    print(output_label + ' laplace')
    print(np.sum(np.array(p_laplace_vec)<0.05)/len(p_laplace_vec))
    print('')
    plt.title(output_label+' pval W uniform wrt laplace')
    plt.show()
    plt.plot(tyrsvec/1e6,p_xJ_vec)
    plt.axhline(y=0.05,color='red')
    print(output_label + ' xJ')
    print(np.sum(np.array(p_xJ_vec)<0.05)/len(p_xJ_vec))
    print('')
    plt.title(output_label+' pval W uniform wrt xJ')
    plt.show()
    plt.plot(tyrsvec/1e6,p_xcc_vec)
    plt.axhline(y=0.05,color='red')
    print(output_label + ' xcc')
    print(np.sum(np.array(p_xcc_vec)<0.05)/len(p_xcc_vec))
    print('')
    plt.title(output_label+' pval W uniform wrt xcc')
    plt.show()
    plt.plot(tyrsvec/1e6,p_xI_vec)
    plt.axhline(y=0.05,color='red')
    print(output_label + ' xinvar')
    print(np.sum(np.array(p_xI_vec)<0.05)/len(p_xI_vec))
    print('')
    plt.title(output_label+' pval W uniform wrt xinvar')
    plt.show()