# -*- coding: utf-8 -*-
"""
Created on %25-Feb-2018

@author: %A.Nishanth C00294860
editted on 12-Oct-2020 at 5.03p.m
"""

import csv
import os
import pickle
import numpy as np
"""
This script is written for checking the MOTIF correlation coefficient among ONGO and TSG

RV_2_coefficient: This doesn't consider the same property varies, thus it behave like correlation coeficient
                    it values are inbetween -1 to 1 
"""
#%% then load the numpy array
#sing_mat_temp = np.load('sing_mat_0.npy')



load_ongo_dir = 'C:/Users/nishy/Desktop/Documents/Fall_2020/codes/Level-3_RV_coefficient/results/temp_results_ONGO_all_positive_prop_represent'
load_tsg_dir = 'C:/Users/nishy/Desktop/Documents/Fall_2020/codes/Level-3_RV_coefficient/results/temp_results_TSG_all_positive_prop_represent'
os.chdir('/')
os.chdir(load_ongo_dir)  
RV_coeff_ongo = np.load('RV_coeff_ongo_unique.npy')
os.chdir('/')
os.chdir(load_tsg_dir)  
RV_coeff_tsg = np.load('RV_coeff_tsg_unique.npy')
# create the numpy array to stor RV_2 coefficient
RV_coeff = np.zeros((len(RV_coeff_ongo),len(RV_coeff_tsg)))    

#%% just calculate the RV_2 coefficient for two matrices
for i in range(0,len(RV_coeff_ongo)):
    os.chdir('/')
    os.chdir(load_ongo_dir)
    sing_mat_main_rv = np.load(''.join(['sing_ongo_mat_' + str(i) , '.npy']))
    for j in range(0,len(RV_coeff_tsg)):
        os.chdir('/')
        os.chdir(load_tsg_dir)
        sing_mat_temp_rv = np.load(''.join(['sing_tsg_mat_' + str(j) , '.npy']))
        numerator_rv = np.sum(np.matmul(np.transpose(sing_mat_main_rv), sing_mat_temp_rv))
        denominater_rv = np.sum((np.matmul(np.transpose(sing_mat_main_rv), sing_mat_main_rv))*np.sum(np.matmul(np.transpose(sing_mat_temp_rv), sing_mat_temp_rv)))
        RV_coeff[i,j] = numerator_rv/(denominater_rv**0.5)
#%%

        #%%
os.chdir('/')
os.chdir(load_ongo_dir)  
RV2_coeff_ongo = np.load('RV2_coeff_ongo_unique.npy')
os.chdir('/')
os.chdir(load_tsg_dir)  
RV2_coeff_tsg = np.load('RV2_coeff_tsg_unique.npy')
# create the numpy array to stor RV_2 coefficient
RV2_coeff = np.zeros((len(RV2_coeff_ongo),len(RV2_coeff_tsg)))    

#%% just calculate the RV_2 coefficient for two matrices
for i in range(0,len(RV2_coeff_ongo)):
    os.chdir('/')
    os.chdir(load_ongo_dir)
    sing_mat_main_rv2 = np.load(''.join(['rv_2_sing_ongo_mat_' + str(i) , '.npy']))
    for j in range(0,len(RV2_coeff_tsg)):
        os.chdir('/')
        os.chdir(load_tsg_dir)
        sing_mat_temp_rv2 = np.load(''.join(['rv_2_sing_tsg_mat_' + str(j) , '.npy']))
        numerator_rv2 = np.sum(np.matmul(np.transpose(sing_mat_main_rv2), sing_mat_temp_rv2))
        denominater_rv2 = np.sum((np.matmul(np.transpose(sing_mat_main_rv2), sing_mat_main_rv2))*np.sum(np.matmul(np.transpose(sing_mat_temp_rv2), sing_mat_temp_rv2)))
        RV2_coeff[i,j] = numerator_rv2/(denominater_rv2**0.5)
#%%
os.chdir('/')
os.chdir('C:/Users/nishy/Desktop/Documents/Fall_2020/codes/Level-3_RV_coefficient/results')
np.save('RV2_coeff_ongo_vs_tsg_all_positive_prop',RV2_coeff)
print("RV coeeficient minimum: ",round(np.min(RV_coeff),3))
print("RV2 coeeficient minimum: ",round(np.min(RV2_coeff),3))