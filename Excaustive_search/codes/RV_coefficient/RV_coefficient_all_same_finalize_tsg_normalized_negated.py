# -*- coding: utf-8 -*-
"""
Created on %25-Feb-2018

@author: %A.Nishanth C00294860
editted on 12-Oct-2020 at 5.38p.m

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
saving_dir = 'C:/Users/nishy/Desktop/Documents/Fall_2020/codes/Level-3_RV_coefficient/results/temp_results_TSG_all_negative_prop_represent'
loading_dir = 'C:/Users/nishy/Desktop/Documents/Fall_2020/codes/Level-3_excaustive_search/Results_and_summary'
os.chdir(loading_dir)
print(os.getcwd())
# then load the text file and get the details of MOTIF
name = "summary_tsg_unique.txt"

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
    '''change made here
    '''
    temp_pass=0
    chk=0
    # then using the tab number to findout the center atom
    while temp_pass<6:
        if temp[chk] == '\t':
            temp_pass=temp_pass+1
        chk=chk+1
    #then group the amino acids from the temp 
    motif_group_temp_1 =[]
    motif_group_temp=[temp[chk]]
    for i in range(chk,len(temp)-1):
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
        if (motif_group_temp[i]  =="M"):
            #        print(motif_group_temp[i]  )
            property_main[i,:]=np.array([-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,0.207070707,0.079276773])
        elif(motif_group_temp[i]  =="R"):
            property_main[i,:]=np.array([1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,1,0.146464646,0.065368567])
        elif(motif_group_temp[i]  =="K"):
            property_main[i,:]=np.array([1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,1,0.747474747,1])
        elif(motif_group_temp[i]  =="D"):
            property_main[i,:]=np.array([1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,-1,-1,0.404040404,0.02364395])
        elif(motif_group_temp[i]  =="E"):
            property_main[i,:]=np.array([1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,-1,-1,0.439393939,0.066759388])
        elif(motif_group_temp[i]  =="Q"):
            property_main[i,:]=np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0.166666667,0.063977747])
        elif(motif_group_temp[i]  =="N"):
            property_main[i,:]=np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0,0.043115438])
        elif(motif_group_temp[i]  =="H"):
            property_main[i,:]=np.array([-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,1,0.085858586,0.009735744])
        elif(motif_group_temp[i]  =="S"):
            property_main[i,:]=np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0.176767677,0.069541029])
        elif(motif_group_temp[i]  =="T"):
            property_main[i,:]=np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0.161616162,0.061196106])
        elif(motif_group_temp[i]  =="Y"):
            property_main[i,:]=np.array([-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,0.156565657,0.068150209])
        elif(motif_group_temp[i]  =="C"):
            property_main[i,:]=np.array([-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,0])
        elif(motif_group_temp[i]  =="W"):
            property_main[i,:]=np.array([-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,0.297979798,0.093184979])
        elif(motif_group_temp[i]  =="A"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.54040404,0.089012517])
        elif(motif_group_temp[i]  =="I"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.484848485,0.084840056])
        elif(motif_group_temp[i]  =="L"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.404040404,0.090403338])
        elif(motif_group_temp[i]  =="F"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,0.222222222,0.121001391])
        elif(motif_group_temp[i]  =="V"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.464646465,0.080667594])
        elif(motif_group_temp[i]  =="P"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0.909090909,0.038942976])
        elif(motif_group_temp[i]  =="G"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0.404040404,0.087621697])
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
    np.save(''.join(['sing_tsg_mat_',str(a-2)]), sing_mat)
    np.save(''.join(['rv_2_sing_tsg_mat_',str(a-2)]),  rv_2_sing)
#%% then load the numpy array
#sing_mat_temp = np.load('sing_mat_0.npy')

# create the numpy array to stor RV_2 coefficient
RV_coeff = np.zeros((len(data)-4,len(data)-4))   
RV2_coeff = np.zeros((len(data)-4,len(data)-4))    

#%% just calculate the RV_2 coefficient for two matrices
for i in range(0,len(data)-4):
    sing_mat_main_rv2 = np.load(''.join(['rv_2_sing_tsg_mat_' + str(i) , '.npy']))
    sing_mat_main = np.load(''.join(['sing_tsg_mat_' + str(i) , '.npy']))
    for j in range(0,len(data)-4):
        sing_mat_temp_rv2 = np.load(''.join(['rv_2_sing_tsg_mat_' + str(j) , '.npy']))
        sing_mat_temp = np.load(''.join(['sing_tsg_mat_' + str(j) , '.npy']))
        numerator = np.sum(np.square(np.matmul(np.transpose(sing_mat_main), sing_mat_temp)))
        denominater = np.sum(np.square(np.matmul(np.transpose(sing_mat_main), sing_mat_main)))*np.sum(np.square(np.matmul(np.transpose(sing_mat_temp), sing_mat_temp)))
        RV_coeff[i,j] = numerator/(denominater**0.5)
        
        numerator_rv2 = np.sum(np.matmul(np.transpose(sing_mat_main_rv2), sing_mat_temp_rv2))
        denominater_rv2 = np.sum((np.matmul(np.transpose(sing_mat_main_rv2), sing_mat_main_rv2))*np.sum(np.matmul(np.transpose(sing_mat_temp_rv2), sing_mat_temp_rv2)))
        RV2_coeff[i,j] = numerator_rv2/(denominater_rv2**0.5)
#%%
np.save('RV2_coeff_tsg_unique',RV2_coeff)
np.save('RV_coeff_tsg_unique',RV_coeff)
print("RV coeeficient minimum: ",round(np.min(RV_coeff),3))
print("RV2 coeeficient minimum: ",round(np.min(RV2_coeff),3))
