#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 12:51:22 2018

@author: c00294860

This script is basically used to find out the which pdb ids fell into the groups at leat once,
and the PDB ids missing fdrom the groups
"""


import os
import pickle
#hashmap_edge=np.zeros((20,20))

#%% load the files with satisfied mot if groups
os.chdir('/')
saving_dir = "C:/Users/nishy/Documents/Projects_UL/mot_if_STUDY/onco_started/finalised_pikle_results_groups/unique_groups"  
os.chdir(saving_dir )
print(os.getcwd())

index = 0

#%% first load the pikle files group details
pikle_files=[]
for l in os.listdir():
    if l.endswith('_total_hits.p'):
        pikle_files.append(l)
#''.join(["amino_acid_",str(a),"_group_",str(m),"_",str(Total_count_group),"_total_hits",".p"]), "wb")) 
#%% then arrage them according to the highest hit to the lowest hit order
# first read the names and take the hit points
hit_points = []
#center_aminoacid_index=[]
for i in range(0, len(pikle_files)):
    name_list = pikle_files[i].split("_")
    hit_points.append(int(name_list[5]))
#    center_aminoacid_index.append(int(name_list[2]))
#%% order them in ascending order index keys and order them accordingly
assending_order_of_group_hit_points_key=sorted(range(len(hit_points)), key=lambda k: hit_points[k])    
# then stack the pikle files decending order
pikle_files_sorted=[]
#center_aminoacid_sorted = []
for i in range(0,len(assending_order_of_group_hit_points_key)):
    pikle_files_sorted.append(pikle_files[assending_order_of_group_hit_points_key[len(assending_order_of_group_hit_points_key)-i-1]])
#    center_aminoacid_sorted.append(center_aminoacid_index[assending_order_of_group_hit_points_key[len(assending_order_of_group_hit_points_key)-i-1]])
#%% then make the summary from the sorted list
PDB_ids_satisfied = []    
for a in range(0,len(pikle_files_sorted)):    
#a = 0
    os.chdir('/')
    os.chdir(saving_dir)
    PDB_set_detail_temp = pickle.load(open(pikle_files_sorted[a], "rb"))
    # to get the which pdb ids satisfy                         
    
    """
    Given the value it have to findout the vertices and center point from that for that group
    from the grop_string name center amino acid name found
    
    The vertices are get from the indexes trace back from the hash map
    
    from the beginning of the indices to 3 indices before containing details about the group 
    """
        
    """
    inorder to get the PDB ids satisfied the conditions
    
    """
    satisfy_group_temp=[] 
    satisfy_group_PDB_ids=[]
    group_sat_PDBs_list_temp = PDB_set_detail_temp[len(PDB_set_detail_temp)-1]
    #% find the . seperation in PDB files and find the index of that first
    # then seperate them to findout the ids
    index_dot=[]
    list_count =0 
    for i in range(0,len(group_sat_PDBs_list_temp)):
        if group_sat_PDBs_list_temp[i]=='.':
            index_dot.append(i)
  
    #%    #Create the summary list from the details build up to now

    os.chdir('/')
    os.chdir("C:/Users/nishy/Desktop/chk_python_summarry")
    
    #%
    name_list = pikle_files_sorted[a].split("_")
    print('create python file part ',index)
    for i in index_dot:
        PDB_ids_satisfied.append(''.join(group_sat_PDBs_list_temp[i+2:i+6]))
    PDB_ids_satisfied = list(set(PDB_ids_satisfied))
    index = index + 1 
#%% creating an text file with the satisfied PDB_ids
f= open(''.join(['summary_ongo_group_satisfied_PDB_ids.txt']),"w+")
#        f.write('Index'+'\t'+'Hits'+'\t'+ 'Center'+ '\t'+ 'Group'+'\r')
for ids in  PDB_ids_satisfied:
    #insert the group details
    f.write(ids+',\r')

print("Done")
f.close() 



