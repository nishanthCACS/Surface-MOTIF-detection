# -*- coding: utf-8 -*-
"""
Created on %%(29-June-2020) at 2.16p.m

@author: %A.Nishanth C00294860
"""

import os
import pickle
from copy import deepcopy
import numpy as np

loading_dir="C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles"
        
os.chdir('/')
os.chdir(loading_dir)    
selected_ONGO_PDBs=pickle.load(open("ONGO_selected_step_2_all_SITE.p", "rb"))  
selected_TSG_PDBs=pickle.load(open("TSG_selected_step_2_all_SITE.p", "rb"))             

#%%%
'''
Make stats depends on the number of residues presented in the MOTIF group
Thus,
Hshtable used to get the indexing of the arrays 
"R","K","D","E","Q","N","H","S","T","Y","C","M","W","A","I","L","F","V","P","G","U"   
one-residue: stats based on their frequency of the residue occurance
    Implementatuion: this can be done by just an array-of 21 aminoacids and their count
Two-residues: stats based on their frequency of the both of them occured together
    Implementatuion: this can be done by just an array-of 21 x 21 array and update(to avoid delays the 21 x 21 is prefered;
    Eventhough hthe  (21 x 21)/2 where diagonal can be there to increased if the residue occured more(memory efficient)

Three residues or higher can be considred as different ways like consider the groups as it is 
Likeget the stat of center atom and their group of amino acid together
way-1:From Hash table get the index (center residue as key); from the hash table get the indexes of the rest of the residues occurance as well
    thus the array of 21 x 21 conatin the occurance_table

And another array conatin the distace addition and divided by their occurance_table find the average distance
        
So in Way-1 Each MOTIF-groups have seperate hashtables with stats of they occurance count
--------------------------------------------------------------------------------------------------
In way-2
the rest of the residues are with count of 
Three-residues: Here onward center residue and the rest group started

Any of the MOTIF groups(contain >3 residues) can be presented as ceter residue combined with other residue groups
'''

Hash_residue_to_index={}
j=0
for key in ["R","K","D","E","Q","N","H","S","T","Y","C","M","W","A","I","L","F","V","P","G","U"]:
    Hash_residue_to_index[key]=j
    j=j+1
    
pickle.dump(Hash_residue_to_index, open("Step_2_Hash_residue_to_index.p", "wb" )) 
step_2_way_1_one_residue_group=np.zeros((21))
step_2_way_1_two_residue_group=np.zeros((21,21))
step_2_way_1_three_residue_group=np.zeros((21,21))
step_2_way_1_four_residue_group=np.zeros((21,21))
step_2_way_1_five_residue_group=np.zeros((21,21))

#to hod the distnce between them
step_2_way_1_two_distance_residue_group=np.zeros((21,21))
step_2_way_1_three_distance_residue_group=np.zeros((21,21))
step_2_way_1_four_distance_residue_group=np.zeros((21,21))
step_2_way_1_five_distance_residue_group=np.zeros((21,21))
#%%
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/ONGO_MOTIF_pickles_V91'
loading_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
os.chdir('/')
os.chdir(loading_dir)
PDB_MOTIF_overall_info_step_2=pickle.load(open("PDB_MOTIF_overall_info_step_2_3H9R.p", "rb"))    
MOTIF_step_2_cosidered= PDB_MOTIF_overall_info_step_2['MOTIF_step_2_cosidered']
taken_MOTIF_summery = MOTIF_step_2_cosidered[0]
   
dist_center_residue=  taken_MOTIF_summery['dist_center_residue']
taken_MOTIF_group_aminoacids = taken_MOTIF_summery['taken_MOTIF_group_aminoacids']
#%%

def find_center_residue(dist_center_residue):
    for i in range(0,len(dist_center_residue)):
        if dist_center_residue[i]==0:
            return i
        
def helper_way_1_stat(dist_center_residue,step_2_way_1_any_residue_group,step_2_way_1_any_residue_dis_group):
    center_res_index=find_center_residue(dist_center_residue)
    index_1=Hash_residue_to_index[taken_MOTIF_group_aminoacids[center_res_index]]

    for k in range(0,len(dist_center_residue)):
        if k !=center_res_index:
            index_2=Hash_residue_to_index[taken_MOTIF_group_aminoacids[k]]
            step_2_way_1_any_residue_group[index_1][index_2]= step_2_way_1_any_residue_group[index_1][index_2]+1
            step_2_way_1_any_residue_dis_group[index_1][index_2]= step_2_way_1_any_residue_dis_group[index_1][index_2]+dist_center_residue[k]
    return step_2_way_1_any_residue_group,step_2_way_1_any_residue_dis_group


#%%
if len(dist_center_residue)==1:
    index_1=Hash_residue_to_index[taken_MOTIF_group_aminoacids[0]]
    step_2_way_1_one_residue_group[index_1]=step_2_way_1_one_residue_group[index_1]+1
elif len(dist_center_residue)==2:
    # just update the both hits
    idex_1=Hash_residue_to_index[taken_MOTIF_group_aminoacids[0]]
    idex_2=Hash_residue_to_index[taken_MOTIF_group_aminoacids[1]]
    if idex_1<idex_2:
        step_2_way_1_one_residue_group[idex_1][idex_2]=step_2_way_1_one_residue_group[idex_1][idex_2]+1
        if dist_center_residue[0]>0:
            step_2_way_1_two_distance_residue_group[idex_1][idex_2]=dist_center_residue[0]
        else:
            step_2_way_1_two_distance_residue_group[idex_1][idex_2]=dist_center_residue[1]
    else:
        step_2_way_1_one_residue_group[idex_2][idex_1]=step_2_way_1_one_residue_group[idex_2][idex_1]+1
        if dist_center_residue[0]>0:
            step_2_way_1_two_distance_residue_group[idex_1][idex_2]=dist_center_residue[0]
        else:
            step_2_way_1_two_distance_residue_group[idex_1][idex_2]=dist_center_residue[1]
elif len(dist_center_residue)==3:
    step_2_way_1_three_residue_group,step_2_way_1_three_distance_residue_group= helper_way_1_stat(dist_center_residue,step_2_way_1_three_residue_group,step_2_way_1_three_distance_residue_group)
elif len(dist_center_residue)==4:
    step_2_way_1_four_residue_group,step_2_way_1_four_distance_residue_group= helper_way_1_stat(dist_center_residue,step_2_way_1_four_residue_group,step_2_way_1_four_distance_residue_group)
elif len(dist_center_residue)==5:
    step_2_way_1_five_residue_group,step_2_way_1_five_distance_residue_group= helper_way_1_stat(dist_center_residue,step_2_way_1_five_residue_group,step_2_way_1_five_distance_residue_group)
else:
    print("some thing here", PDB)
'''
Assigning seperate hash tables for groups higher than 
 
'''
    