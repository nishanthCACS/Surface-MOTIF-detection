# -*- coding: utf-8 -*-
"""
Created on %03-Nov-2020 at 3.27p.m

@author: %A.Nishanth C00294860
"""
import os
import pickle
from Step_2_summerise_small_MOTIF_info_class import  Step_2_group_samples_MOTIF_info


loading_dir="C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles"
        
os.chdir('/')
os.chdir(loading_dir)    
TSG_SITE_satisfied_all_no_thresh =pickle.load(open("TSG_SITE_satisfied_all_no_thresh.p", "rb"))  
ONGO_SITE_satisfied_all_no_thresh =pickle.load(open("ONGO_SITE_satisfied_all_no_thresh.p", "rb"))  

TSG_SITE_satisfied_all_no_thresh_all=list(sum(TSG_SITE_satisfied_all_no_thresh,[]))
ONGO_SITE_satisfied_all_no_thresh_all=list(sum(ONGO_SITE_satisfied_all_no_thresh,[]))
overlapped = pickle.load(open("overlapped_ONGO_TSG_no_thresh_all.p", "rb"))  

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/ONGO_MOTIF_pickles_V91'
saving_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
for pdb in ONGO_SITE_satisfied_all_no_thresh_all:
    if pdb not in overlapped:
        if pdb not in ["4MDQ","721P","6Q3M","5GNK","2RI7","5HYN"]:#,"6P8Q","3ASK","3ASL","4GNE","4GNF","4GNG","4GU0","4GUR"]:# 721P doesn't have signle surface atom due to the surface condition
            obj=Step_2_group_samples_MOTIF_info(loading_dir_MOTIF_pick, saving_dir, pdb)
            if obj.break_cons:
                sel=obj.PDB_MOTIF_overall_info_step_2
                break
#%%