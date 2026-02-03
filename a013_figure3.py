#%%
def ellipse_points_2(a0,b0,siga,sigb,rhoab,probmass):
    import numpy as np
    from scipy.stats import chi2 as sschi2
    cov = np.array([[siga**2,rhoab*siga*sigb],[rhoab*siga*sigb,sigb**2]])
    lenth = 1001
    th = np.linspace(start=0,stop=2*np.pi,num=lenth,endpoint=True)
    costh = np.cos(th)
    sinth=  np.sin(th)
    chi2val = sschi2.ppf(probmass,2)
    # print(chi2val)
    eigenvalues,eigenvectors = np.linalg.eig(cov)
    eigmin = np.min(eigenvalues)
    eigmax = np.max(eigenvalues)
    eigmaxindex = np.argmax(eigenvalues)
    eigmaxvec = eigenvectors[:,eigmaxindex]
    rotation_angle = np.arctan2(eigmaxvec[1],eigmaxvec[0])
    rotation_matrix = np.array([[np.cos(rotation_angle),-np.sin(rotation_angle)],\
                                [np.sin(rotation_angle),np.cos(rotation_angle)]])
    semimajoraxis = np.sqrt(chi2val*eigmax)
    semiminoraxis = np.sqrt(chi2val*eigmin)
    ellipse_x_vec = costh*semimajoraxis
    ellipse_y_vec = sinth*semiminoraxis
    for i in range(lenth):
        xhere = np.array([ellipse_x_vec[i],ellipse_y_vec[i]])
        xhere2 = np.matmul(rotation_matrix.reshape(2,2),xhere.reshape(2,1))
        xherex = xhere2[0][0]
        xherey = xhere2[1][0]
        ellipse_x_vec[i] = xherex
        ellipse_y_vec[i] = xherey
    ellipse_x_vec = ellipse_x_vec + a0
    ellipse_y_vec = ellipse_y_vec + b0
    return ellipse_x_vec,ellipse_y_vec
#%%
def plot_fitted_circles_covariance_3(ax,xcc,ycc,siga,sigb,rhoab,probmass,edgecolor,linestyle_in,lw):
    ellipse_x,ellipse_y = ellipse_points_2(xcc,ycc,siga,sigb,rhoab,probmass)
    ax.plot(ellipse_x,ellipse_y,color=edgecolor,linestyle=linestyle_in,linewidth=lw,alpha=0.5)
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
tyrsvec = np.loadtxt('b003_tyrsvec_' + yrsstr + '.csv',delimiter=',')
aau_mat_planets = np.loadtxt('b003_planets_aau_'+yrsstr+'.csv',delimiter=',')
xinvar = np.loadtxt('b003_planets_xinvar_'+yrsstr+'.csv',delimiter=',')
yinvar = np.loadtxt('b003_planets_yinvar_'+yrsstr+'.csv',delimiter=',')
it = 0
xI = xinvar[it]
yI = yinvar[it]
idegplanets = np.loadtxt('b003_planets_ideg_'+yrsstr+'.csv',delimiter=',')
Wdegplanets = np.loadtxt('b003_planets_nodedeg_'+yrsstr+'.csv',delimiter=',')
sini = np.sin(np.radians(idegplanets))
cosW = np.cos(np.radians(Wdegplanets))
sinW = np.sin(np.radians(Wdegplanets))
qplanets = sini * cosW
pplanets = sini * sinW
xJ = qplanets[0,it]
yJ = pplanets[0,it]
probmasses = [0.68,0.95,0.997]
output_labels = ['BPHSFG','H','S','P']
plot_labels = ['A','H','S','P']
subplotrefs = ['(a) ','(b) ','(c) ','(d) ']
input_conditions = [['Background','Potomac','Hilda','Schubart','Francette','Guinevere'],\
                    ['Hilda'],['Schubart'],['Potomac']]
colors = ['forestgreen','red','blue','goldenrod']
linestyles = ['solid','dotted','dashed','dashdot']
ncons = len(input_conditions)

thetavec = np.linspace(start=0,stop=2*np.pi,num=1001,endpoint=False)
lw_here = 1
alpha_here = 0.8
outfile_plot = 'figure3.png'

fig = plt.figure(figsize=(6,6)) # (width,height)
nrows = 2
ncols = 2
nplots = nrows * ncols
ax = [fig.add_subplot(nrows,ncols,i+1) for i in range(nrows*ncols)]
plt.rcParams.update({'font.size':8})
for iax in range(nrows*ncols):
    a = ax[iax]
    a.tick_params(axis='x',direction='in')
    a.tick_params(axis='y',direction='in')
fig.subplots_adjust(wspace=0.05, hspace=0.05)
# subplot titles go across each row, then on to the next row
subplottitles = ['Hilda group','Hilda family','Schubart family','Potomac family']
for iax in range(nplots):
    ax[iax].text(0.95,0.05,subplotrefs[iax]+subplottitles[iax],horizontalalignment='right',verticalalignment='bottom',\
               transform=ax[iax].transAxes,bbox=dict(facecolor='white',alpha=1))
    ax[iax].grid('on')
for iax in [2,3]:
    ax[iax].set_xlabel(r"$q=\sin(i)\,\cos(\Omega)$")
for iax in [0,2]:
    ax[iax].set_ylabel(r"$p=\sin(i)\,\sin(\Omega)$")
it = 0

xspan = 0.018
xmin = -0.012
ymin = 0.01
xlim = [xmin,xmin+xspan]
ylim = [ymin,ymin+xspan]
xticks = [-0.01,-0.005,0,0.005]
yticks = [0.01,0.015,0.02,0.025]
xticklabels = ['-0.010','-0.005','0','']

for iax in range(nplots):
    for iax2 in range(nplots):
        icon = iax2
        output_label = output_labels[icon]
        plot_label = plot_labels[icon]
        dfvmf = pd.read_csv('b008_'+output_label+'_statistics_vmf_'+yrsstr+'.csv')
        xcc = dfvmf['xcc'][it]
        ycc = dfvmf['ycc'][it]
        ax[iax].scatter(xcc,ycc,color=colors[icon],s=12)
        ax[iax].text(xcc-0.0005,ycc-0.0015,plot_label,color=colors[icon])
    icon = iax
    output_label = output_labels[icon]
    plot_label = plot_labels[icon]
    dfvmf = pd.read_csv('b008_'+output_label+'_statistics_vmf_'+yrsstr+'.csv')
    xcc = dfvmf['xcc'][it]
    ycc = dfvmf['ycc'][it]
    ax[iax].scatter(xcc,ycc,color=colors[icon],s=12)
    angle68rad = np.radians(dfvmf['angle_68_deg'][it])
    angle95rad = np.radians(dfvmf['angle_95_deg'][it])
    angle997rad = np.radians(dfvmf['angle_997_deg'][it])
    dfl_icon = pd.read_csv('b007_laplaceplane_'+output_label+'_'+yrsstr+'.csv')
    dflx = dfl_icon['laplacex'][it]
    dfly = dfl_icon['laplacey'][it]
    ax[iax].scatter(dflx,dfly,color='tomato',linestyle='solid',s=12)
    ax[iax].scatter(xJ,yJ,color='black',s=12)
    ax[iax].scatter(xI,yI,color='brown',s=12)
    ax[iax].text(xJ-0.0005,yJ+0.001,'J',color='black')
    ax[iax].text(xI,yI-0.0015,'I',color='brown')
    ax[iax].text(dflx+0.00025,dfly+0.001,'L',color='tomato')
    x68vec = xcc + np.sin(angle68rad)*np.cos(thetavec)
    y68vec = ycc + np.sin(angle68rad)*np.sin(thetavec)
    ax[iax].plot(x68vec,y68vec,color=colors[icon],linestyle=linestyles[icon],lw=1.5,alpha=0.8)
    x95vec = xcc + np.sin(angle95rad)*np.cos(thetavec)
    y95vec = ycc + np.sin(angle95rad)*np.sin(thetavec)
    ax[iax].plot(x95vec,y95vec,color=colors[icon],linestyle=linestyles[icon],lw=1.5,alpha=0.8)
    x997vec = xcc + np.sin(angle997rad)*np.cos(thetavec)
    y997vec = ycc + np.sin(angle997rad)*np.sin(thetavec)
    ax[iax].plot(x997vec,y997vec,color=colors[icon],linestyle=linestyles[icon],lw=1.5,alpha=0.8)
    ax[iax].set_xlim(xlim)
    ax[iax].set_ylim(ylim)
    ax[iax].set_xticks(xticks)
    ax[iax].set_yticks(yticks)
    ax[iax].set_xticklabels(xticklabels)
    ax[iax].set_box_aspect(1)
    
for iax in [0,1]:
    ax[iax].set_xticklabels([])
for iax in [1,3]:
    ax[iax].set_yticklabels([])

plt.tight_layout()
plt.savefig(outfile_plot,dpi=400)