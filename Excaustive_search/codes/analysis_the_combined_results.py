# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 10:37:36 2017

@author: nishanth
"""

import os
import numpy as np
import pickle
#%%combine the fragmnents of numpy and analysis
os.chdir('/')
#old_chk='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/onco_3.5/pikle_results'
old_chk='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/TSG_3.5/'

os.chdir(old_chk)
unique_all_ids= pickle.load(open("unique_all_ids.p", "rb"))#load the unique ids for analysis
motif_sat_basic_count=np.load('motif_sat_basic_count_all.npy')
motif_sat_basic=np.load('motif_sat_basic_all.npy')
motif_sat_fully=np.load('motif_sat_fully_all.npy')
#chk=pickle.load(open("3vfw.p", "rb"))
#%% analysis the results of combined
# add the row wise data and chk howmany satisy the unique combination
def produce_resuts(dataset):
      print('Mean              :', np.mean(dataset))
      print('Standard deviation:', np.std(dataset))
      print('')
      
motif_sat_basic_added=np.sum(motif_sat_basic,axis= 1)
motif_sat_fully_added=np.sum(motif_sat_fully,axis=1)
motif_sat_basic_count_added=np.sum(motif_sat_basic_count,axis=1)

print("---results-motif_sat_basic_added------")
produce_resuts(motif_sat_basic_added)

print("---results-motif_sat_basic_count_added-----")
produce_resuts(motif_sat_basic_count_added)
print("Maximum hits: ",max(motif_sat_basic_count_added))
print("---results-motif_sat_fully_added------")

produce_resuts(motif_sat_fully_added)

#%%then choose the basics satisfy more than
bins=list(range(0,125,5))
bins.append(335)
hist_basic_count_added, bin_edges = np.histogram(motif_sat_basic_count_added,bins=bins)
hist_fully_added, bin_edges = np.histogram(motif_sat_fully_added,bins=bins)
hist_basic_added, bin_edges = np.histogram(motif_sat_basic_added,bins=bins)
#%%
from matplotlib import pyplot as plt 
bins=list(range(0,125,5))
#bins.append(335)
plt.hist(motif_sat_basic_added,bins) 
#plt.hist(motif_sat_basic_count_added,bins) 
#plt.title("ONGO possible substructure hits histogram") 
#plt.title("TSG possible substructure hits histogram") 

plt.show()

            
#%%
'''
inorder to present an example in the thesis

'''
def val_amino_acid(val):
    """
    this function tells val to amino acid
    """
    if(val==0):
      amino_acid="G"
    elif (val==1):
      amino_acid="M"
    elif(val==2):
      amino_acid="R"
    elif(val==3):
      amino_acid="K"
    elif(val==4):
      amino_acid="D"
    elif(val==5):
      amino_acid="E"
    elif(val==6):
      amino_acid="Q"
    elif(val==7):
      amino_acid="N"
    elif(val==8):
      amino_acid="H"
    elif(val==9):
      amino_acid="S"
    elif(val==10):
      amino_acid="T"
    elif(val==11):
      amino_acid="Y"
    elif(val==12):
      amino_acid="C"
    elif(val==13):
      amino_acid="W"
    elif(val==14):
      amino_acid="A"
    elif(val==15):
      amino_acid="I"
    elif(val==16):
      amino_acid="L"
    elif(val==17):
      amino_acid="F"
    elif(val==18):
      amino_acid="V"
    elif(val==19):
      amino_acid="P"
    return(amino_acid)  
#%%
import collections

compare = lambda x, y: collections.Counter(x) == collections.Counter(y)

chk=list(set(unique_all_ids[0]))
'''
Lets go through and check where the last 
'''
similar_groups=[]
for i in range(0,len(unique_all_ids)):
#for j in range(0,len(unique_all_ids)):
   if compare(list(set(unique_all_ids[i])),chk):
       similar_groups.append(unique_all_ids[i])
#%%
def group_amino(given_list):
    amino_group=[]
    for i in given_list:
       amino_group.append(val_amino_acid(i))
    return amino_group
 
for group in similar_groups:
    print(group_amino(group))