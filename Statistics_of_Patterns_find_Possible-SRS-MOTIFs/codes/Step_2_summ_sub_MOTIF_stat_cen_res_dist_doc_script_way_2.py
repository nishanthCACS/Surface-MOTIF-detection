# -*- coding: utf-8 -*-
"""
Created on %%(29-June-2020) at 2.16p.m

@author: %A.Nishanth C00294860

editted on 14-July-2020 at 2.15pm
"""

import os
import pickle
#from copy import deepcopy
import numpy as np
from Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc_class import Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc

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
Hash_index_to_residue={}
j=0
for key in ["R","K","D","E","Q","N","H","S","T","Y","C","M","W","A","I","L","F","V","P","G","U"]:
    Hash_residue_to_index[key]=j
    Hash_index_to_residue[str(j)]=key
    j=j+1
 
pickle.dump(Hash_residue_to_index, open("Step_2_Hash_residue_to_index.p", "wb" )) 
pickle.dump(Hash_index_to_residue, open("Step_2_Hash_index_to_residue.p", "wb" )) 
#%%
def intialise_stat_way_2_doc(loading_dir_MOTIF_pick):
    '''
    Inorder to sasve the stats intialisation
    '''
    
    saving_dir = ''.join([loading_dir_MOTIF_pick,"/step_2_stat_cent_way_2_doc"])

    os.chdir('/')
    if not os.path.isdir(saving_dir):
        os.makedirs(saving_dir)

    step_2_way_2_three_residue_group=np.zeros((21,21))  
    os.chdir('/')
    os.chdir(saving_dir)             
    for center_residue in ["R","K","D","E","Q","N","H","S","T","Y","C","M","W","A","I","L","F","V","P","G","U"]:
        pickle.dump(step_2_way_2_three_residue_group,open(''.join([center_residue,'_way_2_five_residue_group.p']), "wb" )) 

#%%
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/ONGO_MOTIF_pickles_V91'
intialise_stat_way_2_doc(loading_dir_MOTIF_pick)
for pdb in selected_ONGO_PDBs:
    obj=Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc(Hash_residue_to_index,loading_dir_MOTIF_pick,pdb,way_1=False)
    del obj
    
#%%
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/TSG_MOTIF_pickles_V91'
intialise_stat_way_2_doc(loading_dir_MOTIF_pick)
for pdb in selected_TSG_PDBs:
    obj=Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc(Hash_residue_to_index,loading_dir_MOTIF_pick,pdb,way_1=False)
    del obj
#%%
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/ONGO_MOTIF_pickles_V91'
intialise_stat_way_2_doc(loading_dir_MOTIF_pick)
for pdb in selected_ONGO_PDBs:
    obj=Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc(Hash_residue_to_index,loading_dir_MOTIF_pick,pdb,way_1=False)
    del obj
#%
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/TSG_MOTIF_pickles_V91'
intialise_stat_way_2_doc(loading_dir_MOTIF_pick)
for pdb in selected_TSG_PDBs:
    obj=Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc(Hash_residue_to_index,loading_dir_MOTIF_pick,pdb,way_1=False)
    del obj
#%% to handle errors
'''
Inorder to find the base probability just using the

find the counter of center atom occurance of three groupsbased on the classes
'''
center_count=np.zeros((21,2))

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/ONGO_MOTIF_pickles_V91'
#loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/ONGO_MOTIF_pickles_V91'
loading_dir_ONGO = ''.join([loading_dir_MOTIF_pick,"/step_2_stat_cent_way_2_doc"])

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/TSG_MOTIF_pickles_V91'
#loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/TSG_MOTIF_pickles_V91'
loading_dir_TSG = ''.join([loading_dir_MOTIF_pick,"/step_2_stat_cent_way_2_doc"])

for center_residue in ["R","K","D","E","Q","N","H","S","T","Y","C","M","W","A","I","L","F","V","P","G","U"]:
    os.chdir('/')
    os.chdir(loading_dir_ONGO)
    center_count[Hash_residue_to_index[center_residue]][0]= np.sum(pickle.load(open(''.join([center_residue,'_way_2_five_residue_group.p']), "rb"))) 
    
    os.chdir('/')
    os.chdir(loading_dir_TSG)
    center_count[Hash_residue_to_index[center_residue]][1]= np.sum(pickle.load(open(''.join([center_residue,'_way_2_five_residue_group.p']), "rb")))
#% 
'''
go through the groups again and if it occured 10% of overall group either ONGO or TSG
Then take those groups and their corresponding count of hits of ONGO and TSG

'''
def choose_groups(chek_ONGO,chek_TSG,center_count,Hash_residue_to_index,Hash_index_to_residue):
    chosen_res=[]
    for i in range(0,21):
        for j in range(0,21):
            sel=False
            if chek_ONGO[i][j]>=0.1*center_count[Hash_residue_to_index[center_residue]][0] and center_count[Hash_residue_to_index[center_residue]][0] >0:
                 sel=True
            elif chek_TSG[i][j]>=0.1*center_count[Hash_residue_to_index[center_residue]][1] and center_count[Hash_residue_to_index[center_residue]][1] >0:
                sel=True

            if sel:
                chosen_res.append([Hash_index_to_residue[str(i)],Hash_index_to_residue[str(j)],chek_ONGO[i][j],chek_TSG[i][j]])
    return chosen_res


chosen_res_dic={}
for center_residue in ["R","K","D","E","Q","N","H","S","T","Y","C","M","W","A","I","L","F","V","P","G","U"]:
    os.chdir('/')
    os.chdir(loading_dir_ONGO)
    chek_ONGO= pickle.load(open(''.join([center_residue,'_way_2_five_residue_group.p']), "rb"))
        
    os.chdir('/')
    os.chdir(loading_dir_TSG)
    chek_TSG= pickle.load(open(''.join([center_residue,'_way_2_five_residue_group.p']), "rb"))
    chosen_res_dic[center_residue]=choose_groups(chek_ONGO,chek_TSG,center_count,Hash_residue_to_index,Hash_index_to_residue)
#%%
chosen_res_dic={}
for center_residue in ["C"]:
    os.chdir('/')
    os.chdir(loading_dir_ONGO)
    chek_ONGO= pickle.load(open(''.join([center_residue,'_way_2_five_residue_group.p']), "rb"))
        
    os.chdir('/')
    os.chdir(loading_dir_TSG)
    chek_TSG= pickle.load(open(''.join([center_residue,'_way_2_five_residue_group.p']), "rb"))
    chosen_res_dic[center_residue]=choose_groups(chek_ONGO,chek_TSG,center_count,Hash_residue_to_index,Hash_index_to_residue)