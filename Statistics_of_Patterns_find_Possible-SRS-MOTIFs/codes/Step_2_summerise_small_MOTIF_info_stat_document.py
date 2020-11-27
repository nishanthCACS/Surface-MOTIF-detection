# -*- coding: utf-8 -*-
"""
Created on %%(28-June-2020) at 6.11p.m

@author: %A.Nishanth C00294860

load the pickles created by Step_2_summerise_small_MOTIF_info_class.py and Step_2_summerise_small_MOTIF_info_script.py and summerise the stats
"""
import os
import pickle
from copy import deepcopy
import numpy as np

loading_dir="C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles"
        
os.chdir('/')
os.chdir(loading_dir)    
#TSG_SITE_satisfied_all_no_thresh =pickle.load(open("TSG_SITE_satisfied_all_no_thresh.p", "rb"))  
#ONGO_SITE_satisfied_all_no_thresh =pickle.load(open("ONGO_SITE_satisfied_all_no_thresh.p", "rb"))  
#overlapped = pickle.load(open("overlapped_ONGO_TSG_no_thresh_all.p", "rb"))  
#TSG_SITE_satisfied_all_no_thresh_all=list(sum(TSG_SITE_satisfied_all_no_thresh,[]))
#ONGO_SITE_satisfied_all_no_thresh_all=list(sum(ONGO_SITE_satisfied_all_no_thresh,[]))
selected_ONGO_PDBs=pickle.load(open("ONGO_selected_step_2_all_SITE.p", "rb"))  
selected_TSG_PDBs=pickle.load(open("TSG_selected_step_2_all_SITE.p", "rb"))             

#%%
# the information of overall MOTIF structures and their residues count according to the index

# group_count_info[0]= # of the 1 residue MOTIF groups
# group_count_info[1]= # of the 2 residue MOTIF groups
#     :
# group_count_info[19]= # of the 20 residue MOTIF groups
# group_count_info[20]= # of the residues MOTIF groups > 20
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/ONGO_MOTIF_pickles_V91'
loading_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
os.chdir('/')
os.chdir(loading_dir)
ONGO_group_count_info=np.zeros((21))
for pdb in selected_ONGO_PDBs:
    group_count_info= pickle.load(open( ''.join(["group_count_info_",pdb,".p"]), "rb"))  
    ONGO_group_count_info=ONGO_group_count_info+group_count_info
    
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/TSG_MOTIF_pickles_V91'
loading_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
os.chdir('/')
os.chdir(loading_dir)
TSG_group_count_info=np.zeros((21))
for pdb in selected_TSG_PDBs:
    group_count_info= pickle.load(open( ''.join(["group_count_info_",pdb,".p"]), "rb"))  
    TSG_group_count_info=TSG_group_count_info+group_count_info
#%%
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/ONGO_MOTIF_pickles_V91'
loading_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
os.chdir('/')
os.chdir(loading_dir)
ONGO_group_count_info_all_c_alpha=np.zeros((21))
for pdb in selected_ONGO_PDBs:
    group_count_info= pickle.load(open( ''.join(["group_count_info_",pdb,".p"]), "rb"))  
    ONGO_group_count_info_all_c_alpha=ONGO_group_count_info_all_c_alpha+group_count_info
    
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/TSG_MOTIF_pickles_V91'
loading_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
os.chdir('/')
os.chdir(loading_dir)
TSG_group_count_info_all_c_alpha=np.zeros((21))
for pdb in selected_TSG_PDBs:
    group_count_info= pickle.load(open( ''.join(["group_count_info_",pdb,".p"]), "rb"))  
    TSG_group_count_info_all_c_alpha=TSG_group_count_info_all_c_alpha+group_count_info