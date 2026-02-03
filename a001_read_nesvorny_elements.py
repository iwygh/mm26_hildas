#%%
from astroquery.jplhorizons import Horizons
import pandas as pd
import numpy as np
#%%
infile = 'a000_nesvorny_astorb_astdys_wise_akari_sloan.dat_WITH_PROPER_ELMTS'
outfile = 'b001_astorb_horizons.csv'
top_lines = 0
file1 = open(infile, 'r')
Lines = file1.readlines()
numbers = []
ids = []
nobj = len(Lines) - top_lines
iobj = top_lines
while iobj < len(Lines):
    line = Lines[iobj]
    number = line[0:7]
    idhere = line[7:26]
    while '_' in idhere:
        iii = idhere.index('_')
        idhere = idhere[0:iii] + ' ' + idhere[iii+1:]
    numbers.append(number)
    ids.append(idhere)
    iobj = iobj + 1
jd = 2460796.500000 # May 1, 2025 00:00:00
numbers_list = []
idheres_list = []
targetname_list = [] # elements
H_mag_list = []
datetime_jd_list = []
e_list = []
a_au_list = []
incl_deg_list = []
Omega_deg_list = []
w_deg_list = []
M_deg_list = []
rejected_list = []
for iobj in range(nobj):
    number = numbers[iobj]
    idhere = ids[iobj]
    horizons = Horizons(id=idhere,location='500@10',epochs=jd,id_type='asteroid_name')
    el = horizons.elements()
    if 'spacecraft' not in el['targetname'][0]:
        try: # if Horizons returns a value for everything
            H_mag_list.append(el['H'][0])
            numbers_list.append(number)
            idheres_list.append(idhere)
            targetname_list.append(el['targetname'][0])
            datetime_jd_list.append(el['datetime_jd'][0])
            e_list.append(el['e'][0])
            a_au_list.append(el['a'][0])
            incl_deg_list.append(el['incl'][0])
            Omega_deg_list.append(el['Omega'][0])
            w_deg_list.append(el['w'][0])
            M_deg_list.append(el['M'][0])
            print(iobj+1,nobj,number,idhere,el['a'][0])
        except: # if Horizons screws up anything, don't add to the lists
            rejected_list.append(idhere)
            print(idhere,'rejected')
dictionary = {'number':numbers_list,\
              'idhere':idheres_list,\
              'targetname':targetname_list,\
              'H_mag':H_mag_list,\
              'datetime_jd':datetime_jd_list,\
              'e':e_list,\
              'a_au':a_au_list,\
              'incl_deg':incl_deg_list,\
              'Omega_deg':Omega_deg_list,\
              'w_deg':w_deg_list,\
              'M_deg':M_deg_list}
df = pd.DataFrame.from_dict(dictionary)
df.to_csv(outfile,index=False)
a_au_array = np.array(a_au_list)
print(nobj,len(numbers_list),len(rejected_list))
print(np.min(a_au_array),np.max(a_au_array))
print(rejected_list)
# 3.802579866723511 4.160777908987771
# ['2010 LE43          ']