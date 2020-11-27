# -*- coding: utf-8 -*-
"""
Created on %%(28-June-2020) at 10.33 A.m
editted on 30-Oct-2020 at 5.52p.m

@author: %A.Nishanth C00294860

Using these SITE patterns to get statistics for finding the threshold (highest
 distance residue of substructure from the center residue for each structure; 
 where the residue near to center of gravity is defined as center residue) and 
constraint(new substructures are chosen based on the center residue; 
the center residue’s/amino-acid’s frequency of occurrence in the center residue in  Data-1,
 is considered for selecting center residue consideration)  
for new substructure identification in the available PDBs. 
"""
import os
import pickle
from Step_2_step_1_fixed_summerise_small_MOTIF_info_class import Step_2_summerise_MOTIF_info
loading_dir="C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles"
        
os.chdir('/')
os.chdir(loading_dir)    
TSG_SITE_satisfied_all_no_thresh =pickle.load(open("TSG_SITE_satisfied_all_no_thresh.p", "rb"))  
ONGO_SITE_satisfied_all_no_thresh =pickle.load(open("ONGO_SITE_satisfied_all_no_thresh.p", "rb"))  

TSG_SITE_satisfied_all_no_thresh_all=list(sum(TSG_SITE_satisfied_all_no_thresh,[]))
ONGO_SITE_satisfied_all_no_thresh_all=list(sum(ONGO_SITE_satisfied_all_no_thresh,[]))
overlapped = pickle.load(open("overlapped_ONGO_TSG_no_thresh_all.p", "rb"))  

#%% Define center residue based on the center of gravity of the MOTF substructure
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/ONGO_MOTIF_pickles_V91'
saving_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
for pdb in ONGO_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        if pdb !="4MDQ" and  pdb !="721P":# 721P doesn't have signle surface atom due to the surface condition
            obj=Step_2_summerise_MOTIF_info(loading_dir_MOTIF_pick, saving_dir, pdb)
            del obj

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/TSG_MOTIF_pickles_V91'
saving_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
for pdb in TSG_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        obj=Step_2_summerise_MOTIF_info(loading_dir_MOTIF_pick, saving_dir, pdb)
        del obj        
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/ONGO_MOTIF_pickles_V91'
saving_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
for pdb in ONGO_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        if pdb !="4MDQ" and  pdb !="721P":# 721P doesn't have signle surface atom due to the surface condition
            obj=Step_2_summerise_MOTIF_info(loading_dir_MOTIF_pick, saving_dir, pdb)
            del obj   
            
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/TSG_MOTIF_pickles_V91'
saving_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
for pdb in TSG_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        obj=Step_2_summerise_MOTIF_info(loading_dir_MOTIF_pick, saving_dir, pdb)
        del obj
        
#%%
loading_dir="C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles"
os.chdir('/')
os.chdir(loading_dir)   
selected_ONGO_PDBs=[]
for pdb in ONGO_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        if pdb !="4MDQ" and  pdb !="721P":# 721P doesn't have signle surface atom due to the surface 
            selected_ONGO_PDBs.append(pdb)
pickle.dump(selected_ONGO_PDBs, open("ONGO_selected_step_2_all_SITE.p", "wb" )) 

selected_TSG_PDBs=[]
for pdb in TSG_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        selected_TSG_PDBs.append(pdb)
pickle.dump(selected_TSG_PDBs, open("TSG_selected_step_2_all_SITE.p", "wb" )) 
        