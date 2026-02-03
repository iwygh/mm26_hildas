import numpy as np
import pandas as pd
import rebound
import time
#%%
famstr = 'astorb'
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
t0 = time.time()
infile = 'b001_' + famstr + '_horizons.csv'
df = pd.read_csv(infile)
nobj = df.shape[0]
aau_mat = np.zeros((nobj,nt))
e_mat = np.zeros((nobj,nt))
ideg_mat = np.zeros((nobj,nt))
wdeg_mat = np.zeros((nobj,nt))
Wdeg_mat = np.zeros((nobj,nt))
Mdeg_mat = np.zeros((nobj,nt))
time_list = []
diff_list = []
starting_aau = []
starting_e = []
starting_ideg = []
starting_wdeg = []
starting_Wdeg = []
starting_Mdeg = []
for iobj in range(nobj):
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
    idhere = df['idhere'][iobj]
    aau = df['a_au'][iobj]
    e = df['e'][iobj]
    ideg = df['incl_deg'][iobj]
    wdeg = df['w_deg'][iobj]
    Wdeg = df['Omega_deg'][iobj]
    Mdeg = df['M_deg'][iobj]
    sim.add(primary=sim.particles[0],m=0,a=aau,e=e,inc=np.radians(ideg),\
           omega=np.radians(wdeg),Omega=np.radians(Wdeg),M=np.radians(Mdeg))
    sim.dt = dt_yrs * 2*np.pi # 0.2 yr (2.5% of Hilda orbit)
    sim.N_active = nplanets + 1 # ie, planets plus sun
    sim.move_to_com()
    for it, t in enumerate(tyrsvec):
        sim.integrate(t*2*np.pi,exact_finish_time=True)
        orbits = sim.orbits(primary=sim.particles[0])
        o = orbits[nplanets]
        aau_mat[iobj,it] = o.a
        e_mat[iobj,it] = o.e
        ideg_mat[iobj,it] = np.degrees(o.inc)
        wdeg_mat[iobj,it] = np.mod(np.degrees(o.omega),360)
        Wdeg_mat[iobj,it] = np.mod(np.degrees(o.Omega),360)
        Mdeg_mat[iobj,it] = np.mod(np.degrees(o.M),360)
        sim.move_to_com()
    aau = aau_mat[iobj,0]
    e = e_mat[iobj,0]
    ideg = ideg_mat[iobj,0]
    wdeg = wdeg_mat[iobj,0]
    Wdeg = Wdeg_mat[iobj,0]
    Mdeg = Mdeg_mat[iobj,0]
    starting_aau.append(aau)
    starting_e.append(e)
    starting_ideg.append(ideg)
    starting_wdeg.append(wdeg)
    starting_Wdeg.append(Wdeg)
    starting_Mdeg.append(Mdeg)
    t1 = time.time()
    time_here = (t1-t0)/60
    time_list.append(time_here)
    if iobj == 0:
        time_diff = 0 
        diff_list.append(time_diff)
    else:
        time_diff = time_list[-1] - time_list[-2]
        diff_list.append(time_diff)
    print(famstr,iobj+1,nobj,idhere,time_here,time_diff) 
    del sim
    sim = None
np.savetxt('b004_' + famstr + '_aau_' + yrsstr + '.csv',aau_mat,delimiter=',')
np.savetxt('b004_' + famstr + '_e_' + yrsstr + '.csv',e_mat,delimiter=',')
np.savetxt('b004_' + famstr + '_ideg_' + yrsstr + '.csv',ideg_mat,delimiter=',')
np.savetxt('b004_' + famstr + '_perideg_' + yrsstr + '.csv',wdeg_mat,delimiter=',')
np.savetxt('b004_' + famstr + '_nodedeg_' + yrsstr + '.csv',Wdeg_mat,delimiter=',')
np.savetxt('b004_' + famstr + '_Mdeg_' + yrsstr + '.csv',Mdeg_mat,delimiter=',')
dictionary = {'a_au':starting_aau,'e':starting_e,'ideg':starting_ideg,\
              'perideg':starting_wdeg,'nodedeg':starting_Wdeg,'Mdeg':starting_Mdeg}
df_starting = pd.DataFrame.from_dict(dictionary)
df_starting.to_csv('b004_' + famstr + '_starting_' + yrsstr + '.csv',index=False)
t1 = time.time()
print((t1-t0)/60)
#%%
from matplotlib import pyplot as plt
plt.plot(np.array(diff_list))