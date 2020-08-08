# -*- coding: utf-8 -*-
"""
Created on %24 Apr 2018

@author: %A.Nishanth C00294860
"""
from motif_group_intialize_pikle_resolution import motif_group_pikle_to_property_pikle
import os
import pickle
#%% Initialization for the pikle files loading 
"""
Script on ongo

"""
working_dir = "saving directory.../onco_started/finalised_pikle_results_groups/resolution_dip_ongo_results"
saving_dir =  "saving directory.../onco_started/finalised_pikle_results_groups/property_results_ongo"
#% then load the pikle results and form the dip set for training
"""
This file go through th epikle files created by 
                            motif_group_intialize_pikle
                            
This file only make the dip exmaples
"""
os.chdir('/')
os.chdir(working_dir)
print(os.getcwd())
chk =  os.listdir()

group_files =[]
for l in os.listdir():
    if "_dip_group_ in l":
        group_files.append(l)
        
# when loading the file first choose the same PDB files onlys
pdb_name=[]
for l in group_files:
    pdb_name.append(l[0:4])
pdb_name_unique = list(set(pdb_name))
#% then load the pikle file

for pdb_id in pdb_name_unique: 
    pdb_obj = motif_group_pikle_to_property_pikle(working_dir, saving_dir, group_files, pdb_name_unique,pdb_id)
    pdb_obj.intialize_sub_group_property_matrix()
    pdb_obj.main_fn_property_from_PDB()       
#%   
#os.chdir('/')
#os.chdir(saving_dir)  
#property_chk = pickle.load(open(''.join(["property_all_pdb_",pdb_id,".p"]), "rb")) 

"""
Script on TSG

"""
working_dir = "saving directory.../TSG/finalised_pikle_results_groups/resolution_dip_tsg_results"
saving_dir =  "saving directory.../TSG/finalised_pikle_results_groups/property_results_tsg"
#% then load the pikle results and form the dip set for training
"""
This file go through th epikle files created by 
                            motif_group_intialize_pikle
                            
This file only make the dip exmaples
"""
os.chdir('/')
os.chdir(working_dir)
print(os.getcwd())
chk =  os.listdir()

group_files =[]
for l in os.listdir():
    if "_dip_group_ in l":
        group_files.append(l)
        
# when loading the file first choose the same PDB files onlys
pdb_name=[]
for l in group_files:
    pdb_name.append(l[0:4])
pdb_name_unique = list(set(pdb_name))
#% then load the pikle file

for pdb_id in pdb_name_unique: 
    pdb_obj = motif_group_pikle_to_property_pikle(working_dir, saving_dir, group_files, pdb_name_unique,pdb_id)
    pdb_obj.intialize_sub_group_property_matrix()
    pdb_obj.main_fn_property_from_PDB()       
##%%   
#os.chdir('/')
#os.chdir(saving_dir)  
#property_chk = pickle.load(open(''.join(["property_all_pdb_",pdb_id,".p"]), "rb")) 
