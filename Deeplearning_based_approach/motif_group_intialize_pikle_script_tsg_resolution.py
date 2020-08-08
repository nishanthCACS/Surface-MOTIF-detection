# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %A.Nishanth C00294860
"""

from motif_group_intialize_pikle_resolution import motif_group_intialize_pikle
import os
import pickle
import numpy as np
import copy#to make deep copy

#%% load the resolution and the pdb ids coreesponding to that list
os.chdir('C:/Users/nishy/Desktop/MOTIF_project/')
id_string = pickle.load(open( "id_tsg.p", "rb" ))
res =pickle.load(open( "res_tsg.p", "rb" ))    

#%% to run the script
working_dir = "C:/Users/nishy/Documents/Projects_UL/mot_if_STUDY/TSG"
saving_dir = "C:/Users/nishy/Documents/Projects_UL/mot_if_STUDY/TSG/finalised_pikle_results_groups/resolution_dip_tsg_results"
os.chdir('/')
os.chdir(working_dir)
print(os.getcwd())
# first check the files in the directory
coordinate_files=[]
for l in os.listdir():
    if l.endswith("_surf_atoms.csv"):
        coordinate_files.append(l)
'''
To selecting the file names
'''
for a in range(0,len(coordinate_files)):
#a=151
    os.chdir('/')
    os.chdir(working_dir)
    name_cor=coordinate_files[a]#to open the coordinate file
    #%check whether all names are same pdb id
    id_name_1=name_cor.split("_")
    pdb_name = id_name_1[0]
    print("---------------------///////--------progress number:  ",a)
    print("")
    print("PDB id in progress: ",id_name_1[0])
    # to find the resolution
    for pdb_i in range(0,len(id_string)):
        if id_name_1[0]==id_string[pdb_i][0:4]:
            resolution = res[pdb_i]
            print("Resolution: ",resolution)
    # creating the object
    pdb_obj = motif_group_intialize_pikle(working_dir,saving_dir,pdb_name,resolution)
    #create the surface by triangle
    pdb_obj.creating_surface_by_triangle()
    pdb_obj.choose_dip_bump_struct()#findout the dip and dumps
    pdb_obj.dip_group_amino_acid_pikle()
    #    pdb_obj.group_amino_acid_pikle()  