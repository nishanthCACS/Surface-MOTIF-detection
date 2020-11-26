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
This script is written for checking th eunique ongo and TSG MOTIF here I don't consider overlapping MOTIFs

This script is basically written for calculating the RV_coefficient_2

For the semidefinite matrix is calculated in the biginning,
This matrix calculate the covariance among the properties in group

RV coefficient: consider the value changed between the same property covariance value between two groups
RV_2_coefficient: This doesn't consider the same property varies, thus it behave like correlation coeficient
                    it values are inbetween -1 to 1 
"""
#%% load the MOTIF file and change them into properties
# inorder to do that load the pikle file earlier created and use that 
os.chdir('/')
saving_dir = 'C:/Users/nishy/Desktop/chk_python_summarry/pikle_results/unique_ongo'
loading_dir = 'C:/Users/nishy/Desktop/chk_python_summarry'#"C:/Users/nishy/Documents/Projects_UL/mot_if_STUDY/onco_started/finalised_pikle_results_groups/unique_groups"
os.chdir(loading_dir)
print(os.getcwd())
# then load the text file and get the details of MOTIF
name = "summary_ongo_unique.txt"

with open(name) as f: 
    data = f.readlines() 
#%%
a = 2 # to read the line
for a in range(2,len(data)-1):
    #Then read the center atom from the text file
    temp = data[a]    
    # then using the tab number to findout the center atom
    
    #then group the amino acids from the temp 
    motif_group_temp_1 =[]
    motif_group_temp=[temp[10]]
    for i in range(11,len(temp)-1):
        if temp[i] != '\t':
    #        print(temp[i])
            motif_group_temp_1.append(temp[i])
    # aprt from the center atom cut the other aminoacid by the factor of two 
    # inorder to do that first count the number of appearence and then cut it by factor of two
    ammino_acids_unique = list(set(motif_group_temp_1))  
    ammino_acids_unique_count = []      
    # then count the number of appearence
    for am in ammino_acids_unique:
        count = 0
        for am_in_temp in motif_group_temp_1:
            if am == am_in_temp:
                count = count + 1
        ammino_acids_unique_count.append(count)
    
    #% then create a list of containing amino acids
    for am in range(0,len(ammino_acids_unique)):
       for j in range(0,int(ammino_acids_unique_count[am]/2)):
          motif_group_temp.append(ammino_acids_unique[am])
    
    #% then add the properties of the motifs by the factor of 0.5 of sorrundings
    # then create the numpy array to store the properties
    # since there are 20 properties for each of them
    property_main = np.zeros((len(motif_group_temp),16))
    for i in range(0,len(motif_group_temp)):
    # and the center atom property_main
        if (motif_group_temp[i]=="M"):
    #        print(motif_group_temp[i])
            property_main[i,:]=np.array([0,1,0,0,1,0,1,0,0,0,0,0,1,0,9.21,2.28])
        elif(motif_group_temp[i]=="R"):
            property_main[i,:]=np.array([1,0,0,0,0,1,0,0,0,0,1,0,0,1,9.09,2.18])
        elif(motif_group_temp[i]=="K"):
            property_main[i,:]=np.array([1,0,0,0,0,1,0,0,0,0,1,0,0,1,10.28,8.9])
        elif(motif_group_temp[i]=="D"):
            property_main[i,:]=np.array([1,0,0,0,0,1,0,0,0,1,0,1,0,0,9.6,1.88])
        elif(motif_group_temp[i]=="E"):
            property_main[i,:]=np.array([1,0,0,0,0,1,0,0,0,1,0,1,0,0,9.67,2.19])
        elif(motif_group_temp[i]=="Q"):
            property_main[i,:]=np.array([0,1,0,0,0,1,1,0,0,0,0,0,1,0,9.13,2.17])
        elif(motif_group_temp[i]=="N"):
            property_main[i,:]=np.array([0,1,0,0,0,1,1,0,0,0,0,0,1,0,8.8,2.02])
        elif(motif_group_temp[i]=="H"):
            property_main[i,:]=np.array([0,1,0,0,1,0,0,0,0,0,1,0,0,1,8.97,1.78])
        elif(motif_group_temp[i]=="S"):
            property_main[i,:]=np.array([0,1,0,0,0,1,1,0,0,0,0,0,1,0,9.15,2.21])
        elif(motif_group_temp[i]=="T"):
            property_main[i,:]=np.array([0,1,0,0,0,1,1,0,0,0,0,0,1,0,9.12,2.15])
        elif(motif_group_temp[i]=="Y"):
            property_main[i,:]=np.array([0,1,0,1,0,0,0,1,0,0,0,0,1,0,9.11,2.2])
        elif(motif_group_temp[i]=="C"):
            property_main[i,:]=np.array([0,1,0,0,1,0,1,0,0,0,0,0,1,0,10.78,1.71])
        elif(motif_group_temp[i]=="W"):
            property_main[i,:]=np.array([0,1,0,1,0,0,0,1,0,0,0,0,1,0,9.39,2.38])
        elif(motif_group_temp[i]=="A"):
            property_main[i,:]=np.array([0,0,1,1,0,0,0,0,1,0,0,0,1,0,9.87,2.35])
        elif(motif_group_temp[i]=="I"):
            property_main[i,:]=np.array([0,0,1,1,0,0,0,0,1,0,0,0,1,0,9.76,2.32])
        elif(motif_group_temp[i]=="L"):
            property_main[i,:]=np.array([0,0,1,1,0,0,0,0,1,0,0,0,1,0,9.6,2.36])
        elif(motif_group_temp[i]=="F"):
            property_main[i,:]=np.array([0,0,1,1,0,0,0,1,0,0,0,0,1,0,9.24,2.58])
        elif(motif_group_temp[i]=="V"):
            property_main[i,:]=np.array([0,0,1,1,0,0,0,0,1,0,0,0,1,0,9.72,2.29])
        elif(motif_group_temp[i]=="P"):
            property_main[i,:]=np.array([0,0,1,1,0,0,0,0,0,0,0,0,1,0,10.6,1.99])
        elif(motif_group_temp[i]=="L"):
            property_main[i,:]=np.array([0,0,1,1,0,0,0,0,1,0,0,0,1,0,9.6,2.36])
        elif(motif_group_temp[i]=="G"):
            property_main[i,:]=np.array([0,0,1,1,0,0,0,0,0,0,0,0,1,0,9.6,2.34])
    #% then calculate the covariance(semidefinite) matrix of that
    sing_mat = np.matmul(np.transpose(property_main),property_main)
    rv_2_sing =  sing_mat - np.diag(np.diag(sing_mat))
    os.chdir('/')
    os.chdir(saving_dir)
    """
    Inorder to calculate the RV_coefficient_2 first calculate the covariace matrix 
    for all save as pikle file
    
    then load one by one and calculate the RV_coefficient to the row and save the 
    RV_coeff_2 for all matrices
    """
    np.save(''.join(['sing_ongo_mat_',str(a-2)]), sing_mat)
    np.save(''.join(['rv_2_sing_ongo_mat_',str(a-2)]),  rv_2_sing)
#%% then load the numpy array
#sing_mat_temp = np.load('sing_mat_0.npy')

# create the numpy array to stor RV_2 coefficient
RV_coeff = np.zeros((len(data)-4,len(data)-4))   
RV2_coeff = np.zeros((len(data)-4,len(data)-4))    

#%% just calculate the RV_2 coefficient for two matrices
for i in range(0,len(data)-4):
    sing_mat_main_rv2 = np.load(''.join(['rv_2_sing_ongo_mat_' + str(i) , '.npy']))
    sing_mat_main = np.load(''.join(['sing_ongo_mat_' + str(i) , '.npy']))
    for j in range(0,len(data)-4):
        sing_mat_temp_rv2 = np.load(''.join(['rv_2_sing_ongo_mat_' + str(j) , '.npy']))
        sing_mat_temp = np.load(''.join(['sing_ongo_mat_' + str(j) , '.npy']))
        numerator = np.sum(np.square(np.matmul(np.transpose(sing_mat_main), sing_mat_temp)))
        denominater = np.sum(np.square(np.matmul(np.transpose(sing_mat_main), sing_mat_main)))*np.sum(np.square(np.matmul(np.transpose(sing_mat_temp), sing_mat_temp)))
        RV_coeff[i,j] = numerator/(denominater**0.5)
        
        numerator_rv2 = np.sum(np.matmul(np.transpose(sing_mat_main_rv2), sing_mat_temp_rv2))
        denominater_rv2 = np.sum((np.matmul(np.transpose(sing_mat_main_rv2), sing_mat_main_rv2))*np.sum(np.matmul(np.transpose(sing_mat_temp_rv2), sing_mat_temp_rv2)))
        RV2_coeff[i,j] = numerator_rv2/(denominater_rv2**0.5)
#%%
np.save('RV2_coeff_ongo_unique',RV2_coeff)
np.save('RV_coeff_ongo_unique',RV_coeff)