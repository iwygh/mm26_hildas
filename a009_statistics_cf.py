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
    return ellipse_x_vec,ellipse_y_vec,semimajoraxis,semiminoraxis
#%%
def delta_xyiW(a0,b0,siga,sigb,rhoab,probmass):
    ellipse_x_vec,ellipse_y_vec,semimajoraxis,semiminoraxis = ellipse_points_2(a0,b0,siga,sigb,rhoab,probmass)
    ellipse_i_rad_vec = np.arcsin(np.sqrt(ellipse_x_vec**2+ellipse_y_vec**2))
    sini = np.sin(ellipse_i_rad_vec)
    ellipse_W_rad_vec = np.arctan2(ellipse_y_vec/sini,ellipse_x_vec/sini)
    ellipse_i_deg_vec = np.degrees(ellipse_i_rad_vec)
    ellipse_W_deg_vec = np.degrees(ellipse_W_rad_vec)
    imin = np.min(ellipse_i_deg_vec)
    imax = np.max(ellipse_i_deg_vec)
    delta_i_deg = (imax-imin)/2
    W_1 = np.min(ellipse_W_deg_vec)
    W_2 = np.max(ellipse_W_deg_vec)
    W_1 = np.mod(W_1,360)
    W_2 = np.mod(W_2,360)
    W_a = np.min([W_1,W_2])
    W_b = np.max([W_1,W_2])
    if 0<W_a<90 and 270<W_b<360:
        W_b2 = -(360-W_b)
        Wmin = W_b2
        Wmax = W_a
    else:
        Wmin = W_a
        Wmax = W_b
    delta_W_deg = (Wmax-Wmin)/2
    xmin = np.min(ellipse_x_vec)
    xmax = np.max(ellipse_x_vec)
    ymin = np.min(ellipse_y_vec)
    ymax = np.max(ellipse_y_vec)
    delta_x_deg = np.degrees(np.arcsin((xmax-xmin)/2))
    delta_y_deg = np.degrees(np.arcsin((ymax-ymin)/2))
    quasiradius = np.sqrt(semimajoraxis*semiminoraxis)
    delta_mean_deg = np.degrees(np.arcsin(quasiradius))
    return delta_x_deg,delta_y_deg,delta_i_deg,delta_W_deg,delta_mean_deg
#%%
import numpy as np
import pandas as pd
# import circle_fit as cf
import a000_circlefitcopy as cf
# from circle_fit import standardLSQ
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
ideg_mat_astorb = np.loadtxt('b004_astorb_ideg_'+yrsstr+'.csv',delimiter=',')
Wdeg_mat_astorb = np.loadtxt('b004_astorb_nodedeg_'+yrsstr+'.csv',delimiter=',')
# 19 rows of matrix are:
    # xcc,ycc,icc,Wcc,sigx,sigy,rhoxy,sx1,sx2,sx3,sy1,sy2,sy3,si1,si2,si3,sW1,sW2,sW3
    # all angles in degrees
    # sx,sy,rhoxy - stdevs for x and y; correlation coefficient
    # si1,si2,si3,sW1,sW2,sW3 - vmf 68,95,99.7% uncertainties for icc,Wcc
seed_asteroids = ['Potomac','Hilda','Schubart','Ismene']
nseeds = len(seed_asteroids)
for iseed in range(nseeds):
    seed = seed_asteroids[iseed]
    ideg_mat = np.loadtxt('b005_'+seed+'_ideg_clones_'+yrsstr+'.csv',delimiter=',')
    Wdeg_mat = np.loadtxt('b005_'+seed+'_nodedeg_clones_'+yrsstr+'.csv',delimiter=',')
    nobj = ideg_mat.shape[0]
    print(iseed+1,nseeds)
    irad_mat = np.radians(ideg_mat)
    Wrad_mat = np.radians(Wdeg_mat)
    xmat = np.sin(irad_mat)*np.cos(Wrad_mat)
    ymat = np.sin(irad_mat)*np.sin(Wrad_mat)
    zmat = np.cos(irad_mat)
    xccvec = []
    yccvec = []
    rrvec = []
    for it in range(nt):
        point_coordinates = [xmat[:,it],ymat[:,it]]
        xcc,ycc,rr,sigsig = cf.standardLSQ(np.transpose(point_coordinates))
        # xcc,ycc,rr,sigsig = standardLSQ(np.transpose(point_coordinates))
        xccvec.append(xcc)
        yccvec.append(ycc)
        rrvec.append(rr)
    xccvec = np.array(xccvec)
    yccvec = np.array(yccvec)
    rrvec = np.array(rrvec)
    icc_rad_vec = np.arcsin(np.sqrt(xccvec**2+yccvec**2))
    sini = np.sin(icc_rad_vec)
    Wcc_rad_vec = np.arctan2(yccvec/sini,xccvec/sini)
    icc_deg_vec = np.degrees(icc_rad_vec)
    Wcc_deg_vec = np.degrees(Wcc_rad_vec)
    uhatmat = np.zeros((nobj,nt))
    vhatmat = np.zeros((nobj,nt))
    Whatmat = np.zeros((nobj,3,nt))
    phihatmat = np.zeros((nobj,nt))
    xhatmat = np.zeros((nobj,nt))
    yhatmat = np.zeros((nobj,nt))
    deltahatmat = np.zeros((nobj,nt))
    etahatmat = np.zeros((nobj,nt))
    for it in range(nt):
        for iobj in range(nobj):
            uhatmat[iobj,it] = (xmat[iobj,it]-xccvec[it])/rrvec[it]
            vhatmat[iobj,it] = (ymat[iobj,it]-yccvec[it])/rrvec[it]
            Whatmat[iobj,0,it] = uhatmat[iobj,it]
            Whatmat[iobj,1,it] = vhatmat[iobj,it]
            Whatmat[iobj,2,it] = 1
            phihatmat[iobj,it] = np.arctan2(ymat[iobj,it]-yccvec[it],xmat[iobj,it]-xccvec[it])
            xhatmat[iobj,it] = xccvec[it] + rrvec[it]*np.cos(phihatmat[iobj,it])
            yhatmat[iobj,it] = xccvec[it] + rrvec[it]*np.sin(phihatmat[iobj,it])
            deltahatmat[iobj,it] = xmat[iobj,it] - xhatmat[iobj,it]
            etahatmat[iobj,it] = ymat[iobj,it] - yhatmat[iobj,it]
    sigmasquaredvec = []
    for it in range(nt):
        fish1 = deltahatmat[:,it]
        fish2 = etahatmat[:,it]
        fish3 = np.concatenate((fish1,fish2))
        sigmasquaredvec.append(np.var(fish3))
    WTWhatmat = np.zeros((3,3,nt))
    covhatmat = np.zeros((3,3,nt))
    sigavec = []
    sigbvec = []
    sigrvec = []
    rhoabvec = []
    for it in range(nt):
        try:
            WTWhatmat[:,:,it] = np.matmul(np.transpose(Whatmat[:,:,it]),Whatmat[:,:,it])
            covhatmat[:,:,it] = sigmasquaredvec[it] * np.linalg.inv(WTWhatmat[:,:,it])
            sigavec.append(np.sqrt(covhatmat[0,0,it]))
            sigbvec.append(np.sqrt(covhatmat[1,1,it]))
            sigrvec.append(np.sqrt(covhatmat[2,2,it]))
            rhoabvec.append(covhatmat[0,1,it]/(np.sqrt(covhatmat[0,0,it])*np.sqrt(covhatmat[1,1,it])))
        except:
            print(seed,it,nt,'circle fit nans')
            sigavec.append(np.nan)
            sigbvec.append(np.nan)
            sigrvec.append(np.nan)
            rhoabvec.append(np.nan)
    sigavec = np.array(sigavec)
    sigbvec = np.array(sigbvec)
    sigrvec = np.array(sigrvec)
    rhoabvec = np.array(rhoabvec)
    sx1vec = [0]
    sx2vec = [0]
    sx3vec = [0]
    sy1vec = [0]
    sy2vec = [0]
    sy3vec = [0]
    si1vec = [0]
    si2vec = [0]
    si3vec = [0]
    sW1vec = [0]
    sW2vec = [0]
    sW3vec = [0]
    ss1vec = [0]
    ss2vec = [0]
    ss3vec = [0]
    for it in range(nt-1):
        if xccvec[it+1]**2+yccvec[it+1]**2 >= 1:
            print(seed,it,nt,'icc domain error')
        a0 = xccvec[it+1]
        b0 = yccvec[it+1]
        siga = sigavec[it+1]
        sigb = sigbvec[it+1]
        rhoab = rhoabvec[it+1]
        p0 = probmasses[0]
        delta_x_deg,delta_y_deg,delta_i_deg,delta_W_deg,delta_mean_deg = delta_xyiW(a0,b0,siga,sigb,rhoab,p0)
        sx1vec.append(delta_x_deg)
        sy1vec.append(delta_y_deg)
        si1vec.append(delta_i_deg)
        sW1vec.append(delta_W_deg)
        ss1vec.append(delta_mean_deg)
        p1 = probmasses[1]
        delta_x_deg,delta_y_deg,delta_i_deg,delta_W_deg,delta_mean_deg = delta_xyiW(a0,b0,siga,sigb,rhoab,p1)
        sx2vec.append(delta_x_deg)
        sy2vec.append(delta_y_deg)
        si2vec.append(delta_i_deg)
        sW2vec.append(delta_W_deg)
        ss2vec.append(delta_mean_deg)
        p2 = probmasses[2]
        delta_x_deg,delta_y_deg,delta_i_deg,delta_W_deg,delta_mean_deg = delta_xyiW(a0,b0,siga,sigb,rhoab,p2)
        sx3vec.append(delta_x_deg)
        sy3vec.append(delta_y_deg)
        si3vec.append(delta_i_deg)
        sW3vec.append(delta_W_deg)
        ss3vec.append(delta_mean_deg)
    sx1vec = np.array(sx1vec)
    sx2vec = np.array(sx2vec)
    sx3vec = np.array(sx3vec)
    sy1vec = np.array(sy1vec)
    sy2vec = np.array(sy2vec)
    sy3vec = np.array(sy3vec)
    si1vec = np.array(si1vec)
    si2vec = np.array(si2vec)
    si3vec = np.array(si3vec)
    sW1vec = np.array(sW1vec)
    sW2vec = np.array(sW2vec)
    sW3vec = np.array(sW3vec)
    ss1vec = np.array(ss1vec)
    ss2vec = np.array(ss2vec)
    ss3vec = np.array(ss3vec)
    # xcc,ycc,icc,Wcc,sigx,sigy,rhoxy,sx1,sx2,sx3,sy1,sy2,sy3,si1,si2,si3,sW1,sW2,sW3
    dictionary = {'xcc':xccvec,\
                  'ycc':yccvec,\
                  'ideg_cc':icc_deg_vec,\
                  'Wdeg_cc':Wcc_deg_vec,\
                  'sx':sigavec,\
                  'sy':sigbvec,\
                  'rhoxy':rhoabvec,\
                  'sigx_68_deg':sx1vec,\
                  'sigx_95_deg':sx2vec,\
                  'sigx_997_deg':sx3vec,\
                  'sigy_68_deg':sy1vec,\
                  'sigy_95_deg':sy2vec,\
                  'sigy_997_deg':sy3vec,\
                  'sigi_68_deg':si1vec,\
                  'sigi_95_deg':si2vec,\
                  'sigi_997_deg':si3vec,\
                  'sigW_68_deg':sW1vec,\
                  'sigW_95_deg':sW2vec,\
                  'sigW_997_deg':sW3vec,\
                  'sigs_68_deg':ss1vec,\
                  'sigs_95_deg':ss2vec,\
                  'sigs_997_deg':ss3vec}
    dfd = pd.DataFrame.from_dict(dictionary)
    dfd.to_csv('b009_'+seed+'_clones_statistics_cf_Hmax16.3_'+yrsstr+'.csv',index=None)






