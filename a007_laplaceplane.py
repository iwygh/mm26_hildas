#%% laplace plane helper function
def b32_1_fun(alpha):
    import numpy as np
    from scipy import integrate
    I = integrate.quad(b32_1_integrand,0,2*np.pi,args=(alpha))
    b32_1_result = I[0]/np.pi
    error = I[1]
    return b32_1_result, error
#%% laplace plane helper function
def b32_1_integrand(psi,alpha):
    import numpy as np
    integrand = np.cos(psi)/(1-2*alpha*np.cos(psi)+alpha**2)**(3/2)
    return integrand
#%%
def get_GMdict():
    GM_sun = 1
    GM_mercury = 1/6023600
    GM_venus = 1/408523.71
    GM_earthmoon = 1/328900.56
    GM_mars = 1/3098708
    GM_jupiter = 1/1047.3486
    GM_saturn = 1/3497.898
    GM_uranus = 1/22902.98
    GM_neptune = 1/19412.24
    GMdict = {'sun':GM_sun,'mercury':GM_mercury,'venus':GM_venus,\
              'earth':GM_earthmoon,'mars':GM_mars,'jupiter':GM_jupiter,\
              'saturn':GM_saturn,'uranus':GM_uranus,'neptune':GM_neptune}
    return GMdict
#%% compute laplace plane
def laplace_plane_heliocentric_outerplanetsonly_morestuff(a_sample,wmat_planets,\
    Wmat_planets,Mmat_planets,imat_planets,emat_planets,amat_planets): # inputs already in au, rad
    import numpy as np
    mu_sun = 1.0
    n_sample = np.sqrt(mu_sun/a_sample**3) # rad/s
    GMdict = get_GMdict()
    m_mercury = GMdict['mercury']
    m_venus = GMdict['venus']
    m_earth = GMdict['earth']
    m_mars = GMdict['mars']
    m_rockies = m_mercury + m_venus + m_earth + m_mars
    m_jupiter = GMdict['jupiter']/(mu_sun + m_rockies)
    m_saturn = GMdict['saturn']/(mu_sun + m_rockies)
    m_uranus = GMdict['uranus']/(mu_sun + m_rockies)
    m_neptune = GMdict['neptune']/(mu_sun + m_rockies)
    mu_planets = mu_sun * np.array([m_jupiter,m_saturn,m_uranus,m_neptune])
    N = len(mu_planets)
    a_planets = amat_planets
    i_planets = imat_planets
    i_planets = np.radians(imat_planets)
    W_planets = Wmat_planets
    W_planets = np.radians(Wmat_planets)
    n_planets = np.sqrt(mu_sun/a_planets**3)
    B = np.zeros([N,N])
    alpha = np.zeros([N,N])
    alphabar = np.zeros([N,N])
    # eq 7.128, 7.129
    for j in range(N):
        for k in range(N):
           aj = a_planets[j]
           ak = a_planets[k]
           if aj > ak:
               alpha[j,k] = ak/aj
               alphabar[j,k] = 1
           else:
               alpha[j,k] = aj/ak
               alphabar[j,k] = aj/ak
    # eq 7.134, 7.135
    for j in range(N):
        for k in range(N):
            b32_1_result,error = b32_1_fun(alpha[j,k])
            B[j,k] = 1/4 * mu_planets[k]/(mu_sun+mu_planets[j]) * \
                n_planets[j] * alpha[j,k] * alphabar[j,k] * b32_1_result
    for j in range(N):
        B[j,j] = 0
        for k in range(N):
            if k != j:
                b32_1_result,error = b32_1_fun(alpha[j,k])
                B[j,j] = B[j,j] - n_planets[j] * 1/4 * mu_planets[k] / \
                    (mu_sun+mu_planets[j]) * alpha[j,k] * alphabar[j,k] * \
                    b32_1_result
    I_mat = np.zeros([N,N])
    f_list,Ibar = np.linalg.eig(B) # pg 301 below eq 7.138
    q_planets = np.sin(i_planets)*np.cos(W_planets) # eq 7.19
    p_planets = np.sin(i_planets)*np.sin(W_planets) # eq 7.19
    T_cosgamma = np.linalg.solve(Ibar,q_planets) # eq 7.47
    T_singamma = np.linalg.solve(Ibar,p_planets) # eq 7.47
    T = np.sqrt(T_singamma**2+T_cosgamma**2)
    cosgamma = T_cosgamma/T
    singamma = T_singamma/T
    for i in range(N):
        for j in range(N):
            I_mat[j,i] = Ibar[j,i]*T[i] # eq 7.41
    gamma_array = np.mod(np.arctan2(singamma,cosgamma),2*np.pi)
    alpha_list = []
    alphabar_list = []
    # eq 7.128, 7.129
    for i in range(N):
        if a_planets[i] < a_sample:
            alpha_list.append(a_planets[i]/a_sample)
            alphabar_list.append(1)
        else:
            alpha_list.append(a_sample/a_planets[i])
            alphabar_list.append(a_sample/a_planets[i])
    # eq 7.144
    B_list = []
    for i in range(N):
        b32_1_result,error =  b32_1_fun(alpha_list[i])
        B_here = n_sample/4*mu_planets[i]/mu_sun*alpha_list[i]*alphabar_list[i] * \
            b32_1_result;
        B_list.append(B_here)
    # eq 7.143
    B_scalar = -np.sum(B_list);
    # eq 7.76
    mu_list = [] # not mu as in GM, overloaded notation
    for i in range(N):
        mu_here = 0
        for j in range(N):
            mu_here = mu_here + B_list[j]*I_mat[j,i]
        mu_list.append(mu_here)
    q0 = 0
    p0 = 0
    # eq 7.149, 7.150
    for i in range(N):
        qterm = mu_list[i] / (B_scalar-f_list[i])*np.cos(gamma_array[i])
        pterm = mu_list[i] / (B_scalar-f_list[i])*np.sin(gamma_array[i])
        q0 = q0 - qterm
        p0 = p0 - pterm
    i0 = np.arcsin(np.sqrt(q0**2+p0**2))
    W0 = np.arctan2(p0,q0)
    # i0 = np.mod(i0,2*np.pi)
    W0 = np.mod(W0,2*np.pi)
    f1 = f_list[0]
    f2 = f_list[1]
    f3 = f_list[2]
    f4 = f_list[3]
    mu1 = mu_list[0]
    mu2 = mu_list[1]
    mu3 = mu_list[2]
    mu4 = mu_list[3]
    gm1 = gamma_array[0]
    gm2 = gamma_array[1]
    gm3 = gamma_array[2]
    gm4 = gamma_array[3]
    return q0,p0,i0,W0,B_scalar,f1,f2,f3,f4,mu1,mu2,mu3,mu4,gm1,gm2,gm3,gm4
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
Wrad_mat_planets = np.radians(np.loadtxt('b003_planets_nodedeg_'+yrsstr+'.csv',delimiter=','))
wrad_mat_planets = np.radians(np.loadtxt('b003_planets_perideg_'+yrsstr+'.csv',delimiter=','))
Mrad_mat_planets = np.radians(np.loadtxt('b003_planets_Mdeg_'+yrsstr+'.csv',delimiter=','))
irad_mat_planets = np.radians(np.loadtxt('b003_planets_ideg_'+yrsstr+'.csv',delimiter=','))
e_mat_planets = np.loadtxt('b003_planets_e_'+yrsstr+'.csv',delimiter=',')
aau_mat_planets = np.loadtxt('b003_planets_aau_'+yrsstr+'.csv',delimiter=',')
tyrsvec = np.loadtxt('b003_tyrsvec_' + yrsstr + '.csv',delimiter=',')
nt = len(tyrsvec)
dfa = pd.read_csv('b006_astorb_labels_Hmax16.3.csv')
n_astorb = dfa.shape[0]
probmasses = [0.68,0.95,0.997]
Hmax = 16.3
aau_mat = np.loadtxt('b004_astorb_aau_'+yrsstr+'.csv',delimiter=',')
output_labels = ['P','H','S','BPHSFG']
input_conditions = [['Potomac'],['Hilda'],['Schubart'],\
                    ['Background','Potomac','Hilda','Schubart','Francette','Guinevere']]             
ncons = len(input_conditions)
nobj_running_total = 0
for icon in range(ncons):
    output_label = output_labels[icon]
    # input_groups = input_conditions[icon][0]
    input_groups = input_conditions[icon]
    print(output_label,input_groups)
    astorb_indices = []
    for i_astorb in range(n_astorb):
        if (dfa['H_mag'][i_astorb]<=Hmax) and (dfa['family'][i_astorb] in input_groups) :
            astorb_indices.append(i_astorb)
    nobj = len(astorb_indices)
    nobj_running_total = nobj_running_total + nobj
    print(output_label,nobj,nobj_running_total)
    laplacexvec = []
    laplaceyvec = []
    laplaceivec = []
    laplaceWvec = []
    laplaceBvec = []
    f1vec = []
    f2vec = []
    f3vec = []
    f4vec = []
    mu1vec = []
    mu2vec = []
    mu3vec = []
    mu4vec = []
    gm1vec = []
    gm2vec = []
    gm3vec = []
    gm4vec = []
    for it in range(nt):
        time = tyrsvec[it]*2*np.pi
        aau_vec = aau_mat[astorb_indices,it]
        if output_label in ['P','H','S']:
            aau0 = aau_vec[0]
        else:
            aau0 = np.mean(aau_vec)
        q0,p0,i0,W0,B_scalar,f1,f2,f3,f4,mu1,mu2,mu3,mu4,gm1,gm2,gm3,gm4 = \
            laplace_plane_heliocentric_outerplanetsonly_morestuff(\
            aau0,np.degrees(wrad_mat_planets[:,it]),\
            np.degrees(Wrad_mat_planets[:,it]),np.degrees(Mrad_mat_planets[:,it]),\
            np.degrees(irad_mat_planets[:,it]),e_mat_planets[:,it],aau_mat_planets[:,it])
        laplacexvec.append(q0)
        laplaceyvec.append(p0)
        laplaceivec.append(i0)
        laplaceWvec.append(W0)
        laplaceBvec.append(B_scalar)
        f1vec.append(f1)
        f2vec.append(f2)
        f3vec.append(f3)
        f4vec.append(f4)
        mu1vec.append(mu1)
        mu2vec.append(mu2)
        mu3vec.append(mu3)
        mu4vec.append(mu4)
        gm1vec.append(gm1)
        gm2vec.append(gm2)
        gm3vec.append(gm3)
        gm4vec.append(gm4)
    laplacexvec = np.array(laplacexvec)
    laplaceyvec = np.array(laplaceyvec)
    laplaceivec = np.array(laplaceivec)
    laplaceWvec = np.array(laplaceWvec)
    laplaceBvec = np.array(laplaceBvec)
    f1vec = np.array(f1vec)
    f2vec = np.array(f2vec)
    f3vec = np.array(f3vec)
    f4vec = np.array(f4vec)
    mu1vec = np.array(mu1vec)
    mu2vec = np.array(mu2vec)
    mu3vec = np.array(mu3vec)
    mu4vec = np.array(mu4vec)
    gm1vec = np.array(gm1vec)
    gm2vec = np.array(gm2vec)
    gm3vec = np.array(gm3vec)
    gm4vec = np.array(gm4vec)
    dictionary = {'laplacex':laplacexvec,'laplacey':laplaceyvec,'laplaceB':laplaceBvec,\
                  'laplacei':laplaceivec,'laplaceW':laplaceWvec,\
                  'f1':f1vec,'f2':f2vec,'f3':f3vec,'f4':f4vec,\
                  'mu1':mu1vec,'mu2':mu2vec,'mu3':mu3vec,'mu4':mu4vec,\
                  'gm1':gm1vec,'gm2':gm2vec,'gm3':gm3vec,'gm4':gm4vec}
    dfd = pd.DataFrame.from_dict(dictionary)
    dfd.to_csv('b007_laplaceplane_'+output_label+'_'+yrsstr+'.csv',index=None)
    print(icon+1,ncons,output_label,'end')
    
    