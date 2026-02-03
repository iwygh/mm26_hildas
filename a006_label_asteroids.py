#%%
import pandas as pd
# horizons_file = 'b001_astorb_horizons.csv'
horizons_file = 'b001_astorb_horizons_hildas_july2025.csv'
astorb_file = 'a000_nesvorny_astorb_astdys_wise_akari_sloan.dat_WITH_PROPER_ELMTS'
Hilda_file = 'a000_nesvorny_153_Hilda_family.list_090_WO_INTERLOPERS'
Schubart_file = 'a000_nesvorny_1911_Schubart_family.list_060_WO_INTERLOPERS'
Potomac_file = 'a000_nesvorny_1345_Potomac_family.list_140_WO_INTERLOPERS'
Francette_file = 'a000_nesvorny_1212_Francette_family.list_030_WO_INTERLOPERS'
Guinevere_file = 'a000_nesvorny_2483_Guinevere_family.list_040_WO_INTERLOPERS'
# a2008TG106_file = 'a000_nesvorny_269345_2008TG106_family.list_054_WO_INTERLOPERS'
# leave out 2008 TG106 family because it has only 17 members, 16 of which are also in Schubart family - assign remaining member to Background
family_files = [Hilda_file,Schubart_file,Potomac_file,Francette_file,Guinevere_file]
family_top_lines = [12,12,12,11,11]
# family_top_lines = [11,11,11,10,10]
family_names = ['Hilda','Schubart','Potomac','Francette','Guinevere']
ngroups = len(family_names)
idlists_families = []
idlist2 = []
for igroup in range(ngroups):
    infile = family_files[igroup]
    top_lines = family_top_lines[igroup]
    file1 = open(infile, 'r')
    Lines = file1.readlines()
    ids = []
    nobj = len(Lines) - top_lines
    iobj = top_lines
    while iobj < len(Lines):
        line = Lines[iobj]
        idhere = line[7:26]
        while '_' in idhere:
            iii = idhere.index('_')
            idhere = idhere[0:iii] + ' ' + idhere[iii+1:]
            idhere = idhere.lstrip()
            idhere = idhere.rstrip()
        idhere = idhere.lstrip()
        idhere = idhere.rstrip()
        ids.append(idhere)
        iobj = iobj + 1
    idlists_families.append(ids)
dfa = pd.read_csv(horizons_file)
nobj = dfa.shape[0]
idlist_astorb = dfa['idhere'].to_list()
Hlist_astorb = dfa['H_mag'].to_list()
family_list = []
for iobj in range(nobj):
    idhere = idlist_astorb[iobj]
    idhere2 = idhere.lstrip()
    idhere3 = idhere2.rstrip()
    check = 0
    append = 0
    appendages = []
    for igroup in range(ngroups):
        if idhere3 in idlists_families[igroup]:
            family_list.append(family_names[igroup])
            check = 1
            append = append + 1
            appendages.append(family_names[igroup])
    if check == 0:
        family_list.append('Background')
        append = append + 1
        appendages.append('Background')
    if append != 1:
        print(iobj,idhere3,appendages)
for iobj in range(nobj):
    idhere = idlist_astorb[iobj]
    idhere2 = idhere.rstrip()
    idhere3 = idhere2.lstrip()
    idlist2.append(idhere3)
dictionary = {'idhere':idlist2,\
              'H_mag':Hlist_astorb,\
              'family':family_list}
dfo = pd.DataFrame.from_dict(dictionary)
dfo.to_csv('b006_astorb_labels_Hmax16.3.csv',index=False)
#%%
Hmax = 16.3
all_count = 0
potomac_count = 0
schubart_count = 0
hilda_count = 0
francette_count = 0
guinevere_count = 0
for iobj in range(nobj):
    if Hlist_astorb[iobj] <= Hmax:
        all_count = all_count + 1
        if family_list[iobj] == 'Potomac':
            potomac_count = potomac_count + 1
        if family_list[iobj] == 'Schubart':
            schubart_count = schubart_count + 1
        if family_list[iobj] == 'Hilda':
            hilda_count = hilda_count + 1
        if family_list[iobj] == 'Francette':
            francette_count = francette_count + 1
        if family_list[iobj] == 'Guinevere':
            guinevere_count = guinevere_count + 1
print(all_count,hilda_count,schubart_count,potomac_count,francette_count+guinevere_count,\
      all_count-hilda_count-schubart_count-potomac_count-francette_count-guinevere_count)
#%%
all_small_indices = []
potomac_small_indices = []
hilda_small_indices = []
schubart_small_indices = []
all_small_Hmag = []
potomac_small_Hmag = []
hilda_small_Hmag = []
schubart_small_Hmag = []
all_small_ids = []
potomac_small_ids = []
hilda_small_ids = []
schubart_small_ids = []
for iobj in range(nobj):
    if Hlist_astorb[iobj] <= Hmax:
        all_small_indices.append(iobj)
        all_small_Hmag.append(Hlist_astorb[iobj])
        all_small_ids.append(idlist2[iobj])
        if family_list[iobj] == 'Potomac':
            potomac_small_indices.append(iobj)
            potomac_small_Hmag.append(Hlist_astorb[iobj])
            potomac_small_ids.append(idlist2[iobj])
        if family_list[iobj] == 'Schubart':
            schubart_small_indices.append(iobj)
            schubart_small_Hmag.append(Hlist_astorb[iobj])
            schubart_small_ids.append(idlist2[iobj])
        if family_list[iobj] == 'Hilda':
            hilda_small_indices.append(iobj)
            hilda_small_Hmag.append(Hlist_astorb[iobj])
            hilda_small_ids.append(idlist2[iobj])
all_small_nobj = len(all_small_indices)
potomac_small_nobj = len(potomac_small_indices)
schubart_small_nobj = len(schubart_small_indices)
hilda_small_nobj = len(hilda_small_indices)
print(all_small_nobj,hilda_small_nobj,schubart_small_nobj,potomac_small_nobj)
dictionary = {'idhere':all_small_ids,\
              'H_mag':all_small_Hmag}
dfoa = pd.DataFrame.from_dict(dictionary)
dfoa.to_csv('b006_all_Hmax16.3.csv',index=False)
dictionary = {'idhere':potomac_small_ids,\
              'H_mag':potomac_small_Hmag}
dfop = pd.DataFrame.from_dict(dictionary)
dfop.to_csv('b006_potomac_Hmax16.3.csv',index=False)
dictionary = {'idhere':hilda_small_ids,\
              'H_mag':hilda_small_Hmag}
dfoh = pd.DataFrame.from_dict(dictionary)
dfoh.to_csv('b006_hilda_Hmax16.3.csv',index=False)
dictionary = {'idhere':schubart_small_ids,\
              'H_mag':schubart_small_Hmag}
dfos = pd.DataFrame.from_dict(dictionary)
dfos.to_csv('b006_schubart_Hmax16.3.csv',index=False)

