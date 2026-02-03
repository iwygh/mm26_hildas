import numpy as np
import pandas as pd
import rebound
import time
#%%
obj_list = ['Potomac','Hilda','Schubart','Ismene'] # Potomac, Hilda, Schubart, Ismene
horizons_file = 'b001_astorb_horizons.csv'
nclones = 500
da = 0.01
jd = 2460796.5 # May 1, 2025 00:00:00
time_yrs  = int(2e6)
tstep_yrs = int(2e2)
dt_yrs = 0.2
tyrs_str = '2e6yr'
tstepyrs_str = '2e2yr'
dtyrs_str = '0.2yr'
yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
nt = int(time_yrs/tstep_yrs+1)
tyrsvec = np.loadtxt('b003_tyrsvec_'+yrsstr+'.csv',delimiter=',')
GM_sun = 1
GM_mercury = 1/6023600
GM_venus = 1/408523.71
GM_earthmoon = 1/328900.56
GM_mars = 1/3098708
GM_jupiter = 1/1047.3486
GM_saturn = 1/3497.898
GM_uranus = 1/22902.98
GM_neptune = 1/19412.24
GM_sun = 1 + GM_mercury + GM_venus + GM_earthmoon + GM_mars
GM_jupiter = GM_jupiter / GM_sun
GM_saturn = GM_saturn / GM_sun
GM_uranus = GM_uranus / GM_sun
GM_neptune = GM_neptune / GM_sun
GM_list = [GM_jupiter,GM_saturn,GM_uranus,GM_neptune]
ids = ['5','6','7','8']
nplanets = 4
GM_sun = 1
aau_mat_planets = np.loadtxt('b003_planets_aau_'+yrsstr+'.csv',delimiter=',')
e_mat_planets = np.loadtxt('b003_planets_e_'+yrsstr+'.csv',delimiter=',')
ideg_mat_planets = np.loadtxt('b003_planets_ideg_'+yrsstr+'.csv',delimiter=',')
wdeg_mat_planets = np.loadtxt('b003_planets_perideg_'+yrsstr+'.csv',delimiter=',')
Wdeg_mat_planets = np.loadtxt('b003_planets_nodedeg_'+yrsstr+'.csv',delimiter=',')
Mdeg_mat_planets = np.loadtxt('b003_planets_Mdeg_'+yrsstr+'.csv',delimiter=',')
df = pd.read_csv(horizons_file)
idhere_list = df['idhere'].to_list()
idhere_list_2 = []
for idhere in idhere_list:
    idhere_list_2.append(idhere.rstrip())
time_list = []
diff_list = []
t0 = time.time()
for obj_name in obj_list:
    obj_index = idhere_list_2.index(obj_name)
    aau0 = df['a_au'][obj_index]
    e0 = df['e'][obj_index]
    ideg0 = df['incl_deg'][obj_index]
    Wdeg0 = df['Omega_deg'][obj_index]
    wdeg0 = df['w_deg'][obj_index]
    Mdeg0 = df['M_deg'][obj_index]
    aau0_vec = np.random.uniform(low=-da,high=+da,size=nclones)+aau0
    aau0_vec[0] = aau0
    np.savetxt('b005_'+obj_name+'_aau0_clones.csv',aau0_vec,delimiter=',')
    aau_mat = np.zeros((nclones,nt))
    e_mat = np.zeros((nclones,nt))
    ideg_mat = np.zeros((nclones,nt))
    wdeg_mat = np.zeros((nclones,nt))
    Wdeg_mat = np.zeros((nclones,nt))
    Mdeg_mat = np.zeros((nclones,nt))
    for iclone in range(nclones):
        sim = rebound.Simulation()
        sim.add(m = 1,hash = '0')
        sim.integrator = 'whfast'
        for iplanet in range(nplanets):
            GM = GM_list[iplanet]
            aau = aau_mat_planets[iplanet,0]
            e = e_mat_planets[iplanet,0]
            ideg = ideg_mat_planets[iplanet,0]
            wdeg = wdeg_mat_planets[iplanet,0]
            Wdeg = Wdeg_mat_planets[iplanet,0]
            Mdeg = Mdeg_mat_planets[iplanet,0]
            sim.add(primary=sim.particles[0],m=GM,a=aau,e=e,inc=np.radians(ideg),\
                    omega=np.radians(wdeg),Omega=np.radians(Wdeg),M=np.radians(Mdeg))
            sim.move_to_com()
        aau_mat[iclone,0] = aau0_vec[iclone]
        e_mat[iclone,0] = e0
        ideg_mat[iclone,0] = ideg0
        wdeg_mat[iclone,0] = wdeg0
        Wdeg_mat[iclone,0] = Wdeg0
        Mdeg_mat[iclone,0] = Mdeg0
        sim.add(primary=sim.particles[0],m=0,a=aau0_vec[iclone],e=e0,inc=np.radians(ideg0),\
               omega=np.radians(wdeg0),Omega=np.radians(Wdeg0),M=np.radians(Mdeg0))
        sim.dt = dt_yrs * 2*np.pi # 0.2 yr (2.5% of Hilda orbit)
        sim.N_active = nplanets + 1 # ie, planets plus sun
        sim.move_to_com()
        for it, t in enumerate(tyrsvec):
            # if t > 0:
            sim.integrate(t*2*np.pi,exact_finish_time=True)
            orbits = sim.orbits(primary=sim.particles[0])
            o = orbits[nplanets]
            aau_mat[iclone,it] = o.a
            e_mat[iclone,it] = o.e
            ideg_mat[iclone,it] = np.degrees(o.inc)
            wdeg_mat[iclone,it] = np.mod(np.degrees(o.omega),360)
            Wdeg_mat[iclone,it] = np.mod(np.degrees(o.Omega),360)
            Mdeg_mat[iclone,it] = np.mod(np.degrees(o.M),360)
            sim.move_to_com()
        t1 = time.time()
        time_here = (t1-t0)/60
        time_list.append(time_here)
        if iclone == 0:
            time_diff = 0 
            diff_list.append(time_diff)
        else:
            time_diff = time_list[-1] - time_list[-2]
            diff_list.append(time_diff)
        print(obj_name,iclone+1,nclones,time_here,time_diff)
        del sim
        sim = None
    np.savetxt('b005_' + obj_name + '_aau_clones_' + yrsstr + '.csv',aau_mat,delimiter=',')
    np.savetxt('b005_' + obj_name + '_e_clones_' + yrsstr + '.csv',e_mat,delimiter=',')
    np.savetxt('b005_' + obj_name + '_ideg_clones_' + yrsstr + '.csv',ideg_mat,delimiter=',')
    np.savetxt('b005_' + obj_name + '_perideg_clones_' + yrsstr + '.csv',wdeg_mat,delimiter=',')
    np.savetxt('b005_' + obj_name + '_nodedeg_clones_' + yrsstr + '.csv',Wdeg_mat,delimiter=',')
    np.savetxt('b005_' + obj_name + '_Mdeg_clones_' + yrsstr + '.csv',Mdeg_mat,delimiter=',')
t1 = time.time()
print((t1-t0)/60)
from matplotlib import pyplot as plt
plt.plot(np.array(diff_list))