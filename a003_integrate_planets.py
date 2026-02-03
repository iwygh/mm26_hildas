import numpy as np
import rebound
import time
from astroquery.jplhorizons import Horizons
#%%
jd = 2460796.5 # May 1, 2025 00:00:00
time_yrs  = int(2e6)
tstep_yrs = int(2e2)
dt_yrs = 0.2
tyrs_str = '2e6yr'
tstepyrs_str = '2e2yr'
dtyrs_str = '0.2yr'
yrsstr = tyrs_str + '_' + tstepyrs_str + '_' + dtyrs_str + '_jd' + str(jd)
nt = int(time_yrs/tstep_yrs+1)
tyrsvec = np.linspace(0,time_yrs,nt,endpoint=True)
np.savetxt('b003_tyrsvec_' + yrsstr + '.csv',tyrsvec,delimiter=',')
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
aau_mat_planets = np.zeros((nplanets,nt))
e_mat_planets = np.zeros((nplanets,nt))
ideg_mat_planets = np.zeros((nplanets,nt))
wdeg_mat_planets = np.zeros((nplanets,nt))
Wdeg_mat_planets = np.zeros((nplanets,nt))
Mdeg_mat_planets = np.zeros((nplanets,nt))
GM_sun = 1
starting_sim = rebound.Simulation()
starting_sim.add(m = 1,hash = '0')
starting_sim.integrator = 'whfast'
for iobj in range(nplanets):
    idhere = ids[iobj]
    GM = GM_list[iobj]
    horizons = Horizons(id=idhere,location='500@10',epochs=jd)
    el = horizons.elements()
    aau = el['a'][0]
    e = el['e'][0]
    ideg = el['incl'][0]
    wdeg = np.mod(el['w'][0],360)
    Wdeg = np.mod(el['Omega'][0],360)
    Mdeg = np.mod(el['M'][0],360)
    aau_mat_planets[iobj,0] = aau
    e_mat_planets[iobj,0] = e
    ideg_mat_planets[iobj,0] = ideg
    wdeg_mat_planets[iobj,0] = wdeg
    Wdeg_mat_planets[iobj,0] = Wdeg
    Mdeg_mat_planets[iobj,0] = Mdeg
    starting_sim.add(primary=starting_sim.particles[0],m=GM,a=aau,e=e,inc=np.radians(ideg),\
            omega=np.radians(wdeg),Omega=np.radians(Wdeg),M=np.radians(Mdeg))
    starting_sim.move_to_com()
sim = starting_sim
sim.dt = dt_yrs * 2*np.pi # 0.2 yr (2.5% of Hilda orbit)
print('starting sim')
t0 = time.time()
for it, t in enumerate(tyrsvec):
    # if t > 0:
    sim.integrate(t*2*np.pi,exact_finish_time=True)
    orbits = sim.orbits(primary=sim.particles[0])
    for ip in range(nplanets):
        o = orbits[ip]
        aau_mat_planets[ip,it] = o.a
        e_mat_planets[ip,it] = o.e
        ideg_mat_planets[ip,it] = np.degrees(o.inc)
        wdeg_mat_planets[ip,it] = np.mod(np.degrees(o.omega),360)
        Wdeg_mat_planets[ip,it] = np.mod(np.degrees(o.Omega),360)
        Mdeg_mat_planets[ip,it] = np.mod(np.degrees(o.M),360)
    sim.move_to_com()
t1 = time.time()
print((t1-t0)/60)
del sim
del starting_sim
np.savetxt('b003_planets_aau_' + yrsstr + '.csv',aau_mat_planets,delimiter=',')
np.savetxt('b003_planets_e_' + yrsstr + '.csv',e_mat_planets,delimiter=',')
np.savetxt('b003_planets_ideg_' + yrsstr + '.csv',ideg_mat_planets,delimiter=',')
np.savetxt('b003_planets_perideg_' + yrsstr + '.csv',wdeg_mat_planets,delimiter=',')
np.savetxt('b003_planets_nodedeg_' + yrsstr + '.csv',Wdeg_mat_planets,delimiter=',')
np.savetxt('b003_planets_Mdeg_' + yrsstr + '.csv',Mdeg_mat_planets,delimiter=',')
#%%
x_mat_planets = np.zeros((nplanets+1,nt))
y_mat_planets = np.zeros((nplanets+1,nt))
z_mat_planets = np.zeros((nplanets+1,nt))
vx_mat_planets = np.zeros((nplanets+1,nt))
vy_mat_planets = np.zeros((nplanets+1,nt))
vz_mat_planets = np.zeros((nplanets+1,nt))
GM_sun = 1
starting_sim = rebound.Simulation()
starting_sim.add(m = 1,hash = '0')
starting_sim.integrator = 'whfast'
for iobj in range(nplanets):
    idhere = ids[iobj]
    GM = GM_list[iobj]
    horizons = Horizons(id=idhere,location='500@10',epochs=jd)
    el = horizons.elements()
    aau = el['a'][0]
    e = el['e'][0]
    ideg = el['incl'][0]
    wdeg = np.mod(el['w'][0],360)
    Wdeg = np.mod(el['Omega'][0],360)
    Mdeg = np.mod(el['M'][0],360)
    starting_sim.add(primary=starting_sim.particles[0],m=GM,a=aau,e=e,inc=np.radians(ideg),\
            omega=np.radians(wdeg),Omega=np.radians(Wdeg),M=np.radians(Mdeg))
    starting_sim.move_to_com()
sim = starting_sim
it = 0
for iobj in range(nplanets+1):
    x = sim.particles[ip].x
    y = sim.particles[ip].y
    z = sim.particles[ip].z
    vx = sim.particles[ip].vx
    vy = sim.particles[ip].vy
    vz = sim.particles[ip].vz
    x_mat_planets[ip,it] = x
    y_mat_planets[ip,it] = y
    z_mat_planets[ip,it] = z
    vx_mat_planets[ip,it] = vx
    vy_mat_planets[ip,it] = vy
    vz_mat_planets[ip,it] = vz
sim.dt = dt_yrs * 2*np.pi # 0.2 yr (2.5% of Hilda orbit)
print('starting sim')
t0 = time.time()
for it, t in enumerate(tyrsvec):
    # if t > 0:
    sim.integrate(t*2*np.pi,exact_finish_time=True)
    sim.move_to_com()
    for ip in range(nplanets+1):
        x_mat_planets[ip,it] = sim.particles[ip].x
        y_mat_planets[ip,it] = sim.particles[ip].y
        z_mat_planets[ip,it] = sim.particles[ip].z
        vx_mat_planets[ip,it] = sim.particles[ip].vx
        vy_mat_planets[ip,it] = sim.particles[ip].vy
        vz_mat_planets[ip,it] = sim.particles[ip].vz
    sim.move_to_com()
t1 = time.time()
print((t1-t0)/60)
del sim
del starting_sim
#%%
GM_list = [1,GM_jupiter,GM_saturn,GM_uranus,GM_neptune]
qT_vec = np.zeros(nt)
pT_vec = np.zeros(nt)
sT_vec = np.zeros(nt)
for it in range(nt):
    Lvec = np.array([0,0,0])
    for iplanet in range(nplanets+1):
        x = x_mat_planets[iplanet,it]
        y = y_mat_planets[iplanet,it]
        z = z_mat_planets[iplanet,it]
        vx = vx_mat_planets[iplanet,it]
        vy = vy_mat_planets[iplanet,it]
        vz = vz_mat_planets[iplanet,it]
        rterm = np.array([x,y,z])
        vterm = np.array([vx,vy,vz])
        hterm = np.cross(rterm,vterm)
        Lterm = hterm * GM_list[iplanet]
        Lvec = Lvec + Lterm
    hx = Lvec[0]/np.sqrt(Lvec[0]**2+Lvec[1]**2+Lvec[2]**2)
    hy = Lvec[1]/np.sqrt(Lvec[0]**2+Lvec[1]**2+Lvec[2]**2)
    hz = Lvec[2]/np.sqrt(Lvec[0]**2+Lvec[1]**2+Lvec[2]**2)
    qT = -hy
    pT =  hx
    sT =  hz
    qT_vec[it] = qT
    pT_vec[it] = pT
    sT_vec[it] = sT
np.savetxt('b003_planets_xinvar_' + yrsstr + '.csv',qT_vec,delimiter=',')
np.savetxt('b003_planets_yinvar_' + yrsstr + '.csv',pT_vec,delimiter=',')
from matplotlib import pyplot as plt
plt.plot(tyrsvec,qT_vec)
plt.show()
plt.plot(tyrsvec,pT_vec)
plt.show()
plt.plot(tyrsvec,sT_vec)
plt.show()