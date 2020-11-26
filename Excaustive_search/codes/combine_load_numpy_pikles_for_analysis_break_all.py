# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 09:52:54 2017

@author: nishanth
"""
import os
import numpy as np
import pickle
#%%combine the fragmnents of numpy and analysis
os.chdir('/')
os.chdir('N:/Nishanth_data/study in UL/UL projects/Mot_if_study_python/pikle_results')

results_pikle_files=[]
for l in os.listdir():
    if l.endswith(".p"):
        results_pikle_files.append(l)
unique_all_ids= pickle.load(open("unique_all_ids.p", "rb"))
motif_sat_basic=np.zeros((len(unique_all_ids),len(results_pikle_files)))#this one to store the details about the amino acid set same as the other
motif_sat_basic_count=np.zeros((len(unique_all_ids),len(results_pikle_files)))
motif_sat_fully=np.zeros((len(unique_all_ids),len(results_pikle_files)))#this one to store the eaxctly same group
#%% then combine the results
def combine(part_pikle,motif_sat_basic,motif_sat_basic_count,motif_sat_fully,k,adding_num=20000):
    """
    combine the parts pikle to make the wholoe pikle
    : part_pikle- the pikle part used to assign
    : motif_sat_basic
    : motif_sat_basic_count } pikle parts used to combine
    : motif_sat_fully
    : k- starting point in the new pikle
    """
    motif_sat_basic_temp=np.load(''.join(["motif_sat_basic_",part_pikle,".npy"]))
    motif_sat_basic_count_temp=np.load(''.join(["motif_sat_basic_count_",part_pikle,".npy"]))
    motif_sat_fully_temp=np.load(''.join(["motif_sat_fully_",part_pikle,".npy"]))
    motif_sat_basic[k:k+adding_num,:]=motif_sat_basic_temp
    motif_sat_basic_count[k:k+adding_num,:]=motif_sat_basic_count_temp
    motif_sat_fully[k:k+adding_num,:]=motif_sat_fully_temp
    return motif_sat_basic,motif_sat_basic_count,motif_sat_fully
#%% then combine results
motif_sat_basic,motif_sat_basic_count,motif_sat_fully=combine("part_1",motif_sat_basic,motif_sat_basic_count,motif_sat_fully,0)
motif_sat_basic,motif_sat_basic_count,motif_sat_fully=combine("part_2",motif_sat_basic,motif_sat_basic_count,motif_sat_fully,20000)
motif_sat_basic,motif_sat_basic_count,motif_sat_fully=combine("part_3",motif_sat_basic,motif_sat_basic_count,motif_sat_fully,40000)
motif_sat_basic,motif_sat_basic_count,motif_sat_fully=combine("part_4",motif_sat_basic,motif_sat_basic_count,motif_sat_fully,60000)
motif_sat_basic,motif_sat_basic_count,motif_sat_fully=combine("part_5",motif_sat_basic,motif_sat_basic_count,motif_sat_fully,80000)
motif_sat_basic,motif_sat_basic_count,motif_sat_fully=combine("part_6",motif_sat_basic,motif_sat_basic_count,motif_sat_fully,100000,8346)
#%% save the results as numpy
np.save('motif_sat_basic_count_all.npy',motif_sat_basic_count)  
np.save('motif_sat_basic_all.npy',motif_sat_basic)  
np.save('motif_sat_fully_all.npy',motif_sat_fully)  
#%%
