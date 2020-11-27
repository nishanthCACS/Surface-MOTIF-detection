# -*- coding: utf-8 -*-
"""
Created on %%(26-June-2020) at 10.56 A.m

@author: %A.Nishanth C00294860
"""

from Step_1_extract_the_SITE_info_optilt_surface_class import surface_msms_depth_MOTIF__extract_class_quat
import os
import pickle
#%% first extract the MOTIF information
loading_dir="C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles"
        
os.chdir('/')
os.chdir(loading_dir)    
TSG_SITE_satisfied_all_no_thresh =pickle.load(open("TSG_SITE_satisfied_all_no_thresh.p", "rb"))  
ONGO_SITE_satisfied_all_no_thresh =pickle.load(open("ONGO_SITE_satisfied_all_no_thresh.p", "rb"))  
#%%
'''
check the overlapping PDBs in ONGO and TSG
'''
TSG_SITE_satisfied_all_no_thresh_all=list(sum(TSG_SITE_satisfied_all_no_thresh,[]))
ONGO_SITE_satisfied_all_no_thresh_all=list(sum(ONGO_SITE_satisfied_all_no_thresh,[]))
overlapped=[]
for pdb in TSG_SITE_satisfied_all_no_thresh_all:
    if pdb in ONGO_SITE_satisfied_all_no_thresh_all:
        overlapped.append(pdb)
pickle.dump(overlapped, open( "overlapped_ONGO_TSG_no_thresh_all.p", "wb" ) ) 
#%%
'''
First check the PDB detail in which directory
'''
loading_pikle_dir_1='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Spring_2019/results/msms_pikles/2020_SITE_no_threh_V91'
loading_pikle_dir_2='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Spring_2019/results/msms_pikles/2020_SITE_sat'
loading_pikle_dir_3='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Spring_2019/results/msms_pikles/ONGO'
loading_pikle_dir_4='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Spring_2019/results/msms_pikles/TSG'

files_in_dir_1=os.listdir(loading_pikle_dir_1)
files_in_dir_2=os.listdir(loading_pikle_dir_2)
files_in_dir_3=os.listdir(loading_pikle_dir_3)
files_in_dir_4=os.listdir(loading_pikle_dir_4)

#%%
saving_dir='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/ONGO_MOTIF_pickles_V91'
for pdb in ONGO_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        if pdb !="4MDQ":
            if ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_1:      
                loading_pikle_dir=loading_pikle_dir_1
            elif ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_2:      
                loading_pikle_dir=loading_pikle_dir_2
            elif ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_3:     
                loading_pikle_dir=loading_pikle_dir_3
            else:
                print(pdb,"PDB problem ONGO")

            obj = surface_msms_depth_MOTIF__extract_class_quat(loading_pikle_dir, saving_dir, pdb)
            del obj
       
saving_dir='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/TSG_MOTIF_pickles_V91'
for pdb in TSG_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        if ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_1:      
            loading_pikle_dir=loading_pikle_dir_1
        elif ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_2:      
            loading_pikle_dir=loading_pikle_dir_2
        elif ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_4:     
            loading_pikle_dir=loading_pikle_dir_4
        else:
            print(pdb,"PDB problem TSG")
        obj = surface_msms_depth_MOTIF__extract_class_quat(loading_pikle_dir, saving_dir, pdb)
        del obj
        
#%%
saving_dir='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/ONGO_MOTIF_pickles_V91'
for pdb in ONGO_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        if pdb !="4MDQ":
            if ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_1:      
                loading_pikle_dir=loading_pikle_dir_1
            elif ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_2:      
                loading_pikle_dir=loading_pikle_dir_2
            elif ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_3:     
                loading_pikle_dir=loading_pikle_dir_3
            else:
                print(pdb,"PDB problem ONGO")

            obj = surface_msms_depth_MOTIF__extract_class_quat(loading_pikle_dir, saving_dir, pdb,surf_threh_cond=False)
            del obj
#%
saving_dir='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/TSG_MOTIF_pickles_V91'
for pdb in TSG_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        if ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_1:      
            loading_pikle_dir=loading_pikle_dir_1
        elif ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_2:      
            loading_pikle_dir=loading_pikle_dir_2
        elif ''.join(['amino_acid_',pdb,'.p']) in files_in_dir_4:     
            loading_pikle_dir=loading_pikle_dir_4
        else:
            print(pdb,"PDB problem TSG")
        obj = surface_msms_depth_MOTIF__extract_class_quat(loading_pikle_dir, saving_dir, pdb,surf_threh_cond=False)
        del obj