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
                    
Here the edges properties are only considered
thus if the edge contain: ST
then the sum of propererty S and T is considered here
"""
#%% load the MOTIF file and change them into properties
# inorder to do that load the pikle file earlier created and use that 
os.chdir('/')
saving_dir = 'C:/Users/nishy/Desktop/chk_python_summarry/pikle_results/unique_ongo_vertices_sum'
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
    
    #% then create the egdes group inorder to do that crete 
    #two vectors one containing one part of edge vertice
    #      one contain vertice
    #      then another vertice
    edge_1 = []
    edge_2 = []
    # create the center one edges
    for i in range(1,len(motif_group_temp)):
        edge_1.append(motif_group_temp[0])
        edge_2.append(motif_group_temp[i])
    
    #% then add the other vertices
    ver = 1 # this for point which vertice is done and conditioning to select other vertice
    for i in range(13,len(temp)-1):
        if temp[i] != '\t':    
            if ver == 1:
    #            print("ver if ver 1: ", ver)
                edge_1.append(temp[i])
                ver = 2 
    #            print("ver if ver 1: ", ver)
            elif ver == 2:
                edge_2.append(temp[i])
                ver = 1
    #% then add the properties of the motifs by the factor of 0.5 of sorrundings
    # then create the numpy array to store the properties
    # since there are 20 properties for each of them
    property_main = np.zeros((len(edge_1),16))
    for i in range(0,len(edge_1)):
    # and the center atom property_main
        if (edge_1[i]=="M"):
    #        print(edge_1[i])
            property_main[i,:]=np.array([-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,0.207070707,0.079276773])
        elif(edge_1[i]=="R"):
            property_main[i,:]=np.array([1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,1,0.146464646,0.065368567])
        elif(edge_1[i]=="K"):
            property_main[i,:]=np.array([1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,1,0.747474747,1])
        elif(edge_1[i]=="D"):
            property_main[i,:]=np.array([1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,-1,-1,0.404040404,0.02364395])
        elif(edge_1[i]=="E"):
            property_main[i,:]=np.array([1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,-1,-1,0.439393939,0.066759388])
        elif(edge_1[i]=="Q"):
            property_main[i,:]=np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0.166666667,0.063977747])
        elif(edge_1[i]=="N"):
            property_main[i,:]=np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0,0.043115438])
        elif(edge_1[i]=="H"):
            property_main[i,:]=np.array([-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,1,0.085858586,0.009735744])
        elif(edge_1[i]=="S"):
            property_main[i,:]=np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0.176767677,0.069541029])
        elif(edge_1[i]=="T"):
            property_main[i,:]=np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0.161616162,0.061196106])
        elif(edge_1[i]=="Y"):
            property_main[i,:]=np.array([-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,0.156565657,0.068150209])
        elif(edge_1[i]=="C"):
            property_main[i,:]=np.array([-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,0])
        elif(edge_1[i]=="W"):
            property_main[i,:]=np.array([-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,0.297979798,0.093184979])
        elif(edge_1[i]=="A"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.54040404,0.089012517])
        elif(edge_1[i]=="I"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.484848485,0.084840056])
        elif(edge_1[i]=="L"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.404040404,0.090403338])
        elif(edge_1[i]=="F"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,0.222222222,0.121001391])
        elif(edge_1[i]=="V"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.464646465,0.080667594])
        elif(edge_1[i]=="P"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0.909090909,0.038942976])
        elif(edge_1[i]=="G"):
            property_main[i,:]=np.array([-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0.404040404,0.087621697])
    #% then add the property of the other vertice
    for i in range(0,len(edge_2)):
    # and the center atom property_main
        if (edge_1[i]=="M"):
    #        print(edge_1[i])
             property_main[i,:] = property_main[i,:] + np.array([-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,0.207070707,0.079276773])
        elif(edge_1[i]=="R"):
             property_main[i,:] = property_main[i,:] + np.array([1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,1,0.146464646,0.065368567])
        elif(edge_1[i]=="K"):
             property_main[i,:] = property_main[i,:] + np.array([1,-1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,-1,1,0.747474747,1])
        elif(edge_1[i]=="D"):
             property_main[i,:] = property_main[i,:] + np.array([1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,-1,-1,0.404040404,0.02364395])
        elif(edge_1[i]=="E"):
             property_main[i,:] = property_main[i,:] + np.array([1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,-1,-1,0.439393939,0.066759388])
        elif(edge_1[i]=="Q"):
             property_main[i,:] = property_main[i,:] + np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0.166666667,0.063977747])
        elif(edge_1[i]=="N"):
             property_main[i,:] = property_main[i,:] + np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0,0.043115438])
        elif(edge_1[i]=="H"):
             property_main[i,:] = property_main[i,:] + np.array([-1,1,-1,-1,1,-1,-1,-1,-1,-1,1,-1,-1,1,0.085858586,0.009735744])
        elif(edge_1[i]=="S"):
             property_main[i,:] = property_main[i,:] + np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0.176767677,0.069541029])
        elif(edge_1[i]=="T"):
             property_main[i,:] = property_main[i,:] + np.array([-1,1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,1,-1,0.161616162,0.061196106])
        elif(edge_1[i]=="Y"):
             property_main[i,:] = property_main[i,:] + np.array([-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,0.156565657,0.068150209])
        elif(edge_1[i]=="C"):
             property_main[i,:] = property_main[i,:] + np.array([-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,0])
        elif(edge_1[i]=="W"):
             property_main[i,:] = property_main[i,:] + np.array([-1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,0.297979798,0.093184979])
        elif(edge_1[i]=="A"):
             property_main[i,:] = property_main[i,:] + np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.54040404,0.089012517])
        elif(edge_1[i]=="I"):
             property_main[i,:] = property_main[i,:] + np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.484848485,0.084840056])
        elif(edge_1[i]=="L"):
             property_main[i,:] = property_main[i,:] + np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.404040404,0.090403338])
        elif(edge_1[i]=="F"):
             property_main[i,:] = property_main[i,:] + np.array([-1,-1,1,1,-1,-1,-1,1,-1,-1,-1,-1,1,-1,0.222222222,0.121001391])
        elif(edge_1[i]=="V"):
             property_main[i,:] = property_main[i,:] + np.array([-1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,0.464646465,0.080667594])
        elif(edge_1[i]=="P"):
             property_main[i,:] = property_main[i,:] + np.array([-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0.909090909,0.038942976])
        elif(edge_1[i]=="G"):
             property_main[i,:] = property_main[i,:] + np.array([-1,-1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,0.404040404,0.087621697])         
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
    np.save(''.join(['sing_ongo_mat_sum_',str(a-2)]), sing_mat)
    np.save(''.join(['rv_2_sing_ongo_mat_sum_',str(a-2)]),  rv_2_sing)
#%% then load the numpy array
#sing_mat_temp = np.load('sing_mat_0.npy')

# create the numpy array to stor RV_2 coefficient
RV_coeff = np.zeros((len(data)-4,len(data)-4))   
RV2_coeff = np.zeros((len(data)-4,len(data)-4))    

#%% just calculate the RV_2 coefficient for two matrices
for i in range(0,len(data)-4):
    sing_mat_main_rv2 = np.load(''.join(['rv_2_sing_ongo_mat_sum_' + str(i) , '.npy']))
    sing_mat_main = np.load(''.join(['sing_ongo_mat_sum_' + str(i) , '.npy']))
    for j in range(0,len(data)-4):
        sing_mat_temp_rv2 = np.load(''.join(['rv_2_sing_ongo_mat_sum_' + str(j) , '.npy']))
        sing_mat_temp = np.load(''.join(['sing_ongo_mat_sum_' + str(j) , '.npy']))
        numerator = np.sum(np.square(np.matmul(np.transpose(sing_mat_main), sing_mat_temp)))
        denominater = np.sum(np.square(np.matmul(np.transpose(sing_mat_main), sing_mat_main)))*np.sum(np.square(np.matmul(np.transpose(sing_mat_temp), sing_mat_temp)))
        RV_coeff[i,j] = numerator/(denominater**0.5)
        
        numerator_rv2 = np.sum(np.matmul(np.transpose(sing_mat_main_rv2), sing_mat_temp_rv2))
        denominater_rv2 = np.sum((np.matmul(np.transpose(sing_mat_main_rv2), sing_mat_main_rv2))*np.sum(np.matmul(np.transpose(sing_mat_temp_rv2), sing_mat_temp_rv2)))
        RV2_coeff[i,j] = numerator_rv2/(denominater_rv2**0.5)
#%%
np.save('RV2_coeff_sum_ongo_unique',RV2_coeff)
np.save('RV_coeff_sum_ongo_unique',RV_coeff)