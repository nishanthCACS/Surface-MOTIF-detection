# -*- coding: utf-8 -*-
"""
Created on %25-Feb-2018

@author: %A.Nishanth C00294860
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



load_ongo_dir = 'C:/Users/nishy/Desktop/chk_python_summarry/pikle_results/unique_ongo'
load_tsg_dir = 'C:/Users/nishy/Desktop/chk_python_summarry/pikle_results/unique_TSG'

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
os.chdir('C:/Users/nishy/Desktop/chk_python_summarry/pikle_results')
np.save('RV2_coeff_ongo_vs_tsg',RV2_coeff)