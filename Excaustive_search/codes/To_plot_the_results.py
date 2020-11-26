# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %A.Nishanth C00294860
"""

from motif_group_intialize_pikle import motif_group_intialize_pikle
import os
import pickle
import numpy as np
import copy#to make deep copy
#%% to run the script

working_dir = "C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Fall_2018/results/ONGO/test"
saving_dir = "C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Fall_2018/results/ONGO/test"

os.chdir('/')
os.chdir(working_dir)
print(os.getcwd())
# first check the files in the directory
coordinate_files=[]
for l in os.listdir():
    if l.endswith("_surf_atoms.csv"):
        coordinate_files.append(l)
#%%
#a=0 
a=1
os.chdir('/')
os.chdir(working_dir)
name_cor=coordinate_files[a]#to open the coordinate file
#%check whether all names are same pdb id
id_name_1=name_cor.split("_")
pdb_name = id_name_1[0]
print("---------------------///////--------progress number:  ",a)
print("")
print("PDB id in progress: ",id_name_1[0])
# creating the object
pdb_obj = motif_group_intialize_pikle(working_dir,saving_dir,pdb_name)
#create the surface by triangle
pdb_obj.creating_surface_by_triangle()
pdb_obj.group_amino_acid_pikle()
pdb_obj.plot_surface()
