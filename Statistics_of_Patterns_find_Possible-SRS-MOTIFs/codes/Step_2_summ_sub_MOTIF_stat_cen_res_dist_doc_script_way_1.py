# -*- coding: utf-8 -*-
"""
Created on %%(29-June-2020) at 2.16p.m

@author: %A.Nishanth C00294860

updated on 30-Oct-2020
"""

import os
import pickle
from copy import deepcopy
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
j=0
for key in ["R","K","D","E","Q","N","H","S","T","Y","C","M","W","A","I","L","F","V","P","G","U"]:
    Hash_residue_to_index[key]=j
    j=j+1
    
pickle.dump(Hash_residue_to_index, open("Step_2_Hash_residue_to_index.p", "wb" )) 

def intialise_stat_way_1_doc(loading_dir_MOTIF_pick):
    '''
    Inorder to sasve the stats intialisation
    '''
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
    
    
    saving_dir = ''.join([loading_dir_MOTIF_pick,"/step_2_stat_cent_way_1_doc"])

    os.chdir('/')
    if not os.path.isdir(saving_dir):
        os.makedirs(saving_dir)

    os.chdir('/')
    os.chdir(saving_dir)
    pickle.dump(step_2_way_1_one_residue_group, open('step_2_way_1_one_residue_group.p', "wb" )) 
    pickle.dump(step_2_way_1_two_residue_group, open('step_2_way_1_two_residue_group.p', "wb" )) 
    pickle.dump(step_2_way_1_three_residue_group, open('step_2_way_1_three_residue_group.p', "wb" )) 
    pickle.dump(step_2_way_1_four_residue_group, open('step_2_way_1_four_residue_group.p', "wb" )) 
    pickle.dump(step_2_way_1_five_residue_group, open('step_2_way_1_five_residue_group.p', "wb" )) 
    
    pickle.dump(step_2_way_1_two_distance_residue_group, open('step_2_way_1_two_distance_residue_group.p', "wb" )) 
    pickle.dump(step_2_way_1_three_distance_residue_group, open('step_2_way_1_three_distance_residue_group.p', "wb" )) 
    pickle.dump(step_2_way_1_four_distance_residue_group, open('step_2_way_1_four_distance_residue_group.p', "wb" )) 
    pickle.dump(step_2_way_1_five_distance_residue_group, open('step_2_way_1_five_distance_residue_group.p', "wb" )) 

#%%
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/ONGO_MOTIF_pickles_V91'
intialise_stat_way_1_doc(loading_dir_MOTIF_pick)
for pdb in selected_ONGO_PDBs:
    obj=Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc(Hash_residue_to_index,loading_dir_MOTIF_pick,pdb)
    del obj

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/TSG_MOTIF_pickles_V91'
intialise_stat_way_1_doc(loading_dir_MOTIF_pick)
for pdb in selected_TSG_PDBs:
    obj=Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc(Hash_residue_to_index,loading_dir_MOTIF_pick,pdb)
    del obj
#%%
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/ONGO_MOTIF_pickles_V91'
intialise_stat_way_1_doc(loading_dir_MOTIF_pick)
for pdb in selected_ONGO_PDBs:
    obj=Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc(Hash_residue_to_index,loading_dir_MOTIF_pick,pdb)
    del obj
#%
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/TSG_MOTIF_pickles_V91'
intialise_stat_way_1_doc(loading_dir_MOTIF_pick)
for pdb in selected_TSG_PDBs:
    obj=Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc(Hash_residue_to_index,loading_dir_MOTIF_pick,pdb)
    del obj
#%% 
'''to handle errors'''
loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/ONGO_MOTIF_pickles_V91'

loading_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
os.chdir('/')
os.chdir(loading_dir)
PDB_MOTIF_overall_info_step_2 =pickle.load(open(''.join(["PDB_MOTIF_overall_info_step_2_",'4I51',".p"]), "rb"))  
MOTIF_step_2_cosidered= PDB_MOTIF_overall_info_step_2['MOTIF_step_2_cosidered']

for i in range(0,len(MOTIF_step_2_cosidered)):
    taken_MOTIF_summery=MOTIF_step_2_cosidered[i]
#    dist_center_residue=  taken_MOTIF_summery['dist_center_residue']
    taken_MOTIF_group_aminoacids = taken_MOTIF_summery['taken_MOTIF_group_aminoacids']   
    for key in taken_MOTIF_group_aminoacids:
        print(i,'  : ',key)
        print(Hash_residue_to_index[key])
        
#%%
'''
Load the statistics and finalise the results
'''
import pickle
import os

def load_step_2_way_1_results(loading_dir_MOTIF_pick):
    load_dir = ''.join([loading_dir_MOTIF_pick,"/step_2_stat_cent_way_1_doc"])

    os.chdir('/')
    os.chdir(load_dir)
    results={}   
    results["step_2_way_1_one_residue_group"]=pickle.load(open("step_2_way_1_one_residue_group.p", "rb"))  
    results["step_2_way_1_two_residue_group"]=pickle.load(open("step_2_way_1_two_residue_group.p", "rb"))  
    results["step_2_way_1_three_residue_group"]=pickle.load(open("step_2_way_1_three_residue_group.p", "rb"))  
    results["step_2_way_1_four_residue_group"]=pickle.load(open("step_2_way_1_four_residue_group.p", "rb"))  
    results["step_2_way_1_five_residue_group"]=pickle.load(open("step_2_way_1_five_residue_group.p", "rb"))  
        
    results["step_2_way_1_two_distance_residue_group"]=pickle.load(open("step_2_way_1_two_distance_residue_group.p", "rb"))  
    results["step_2_way_1_three_distance_residue_group"]=pickle.load(open("step_2_way_1_three_distance_residue_group.p", "rb"))  
    results["step_2_way_1_four_distance_residue_group"]=pickle.load(open("step_2_way_1_four_distance_residue_group.p", "rb"))  
    results["step_2_way_1_five_distance_residue_group"]=pickle.load(open("step_2_way_1_five_distance_residue_group.p", "rb"))  
    return results

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/ONGO_MOTIF_pickles_V91'
ONGO_MOTIF_soft_results=load_step_2_way_1_results(loading_dir_MOTIF_pick)

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/TSG_MOTIF_pickles_V91'
TSG_MOTIF_soft_results=load_step_2_way_1_results(loading_dir_MOTIF_pick)

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/ONGO_MOTIF_pickles_V91'
ONGO_MOTIF_all_c_alpha_results=load_step_2_way_1_results(loading_dir_MOTIF_pick)

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/TSG_MOTIF_pickles_V91'
TSG_MOTIF_all_c_alpha_results=load_step_2_way_1_results(loading_dir_MOTIF_pick)


#%%
residues_for_index= ["R","K","D","E","Q","N","H","S","T","Y","C","M","W","A","I","L","F","V","P","G","U"]
# take the one rediue
def compare_key_output(ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,compare_key):
    
    differ_soft=ONGO_MOTIF_soft_results[compare_key]-TSG_MOTIF_soft_results[compare_key]
    differ_all_c_alphat=ONGO_MOTIF_all_c_alpha_results[compare_key]-TSG_MOTIF_all_c_alpha_results[compare_key]

    return [ONGO_MOTIF_soft_results[compare_key],TSG_MOTIF_soft_results[compare_key],ONGO_MOTIF_all_c_alpha_results[compare_key],TSG_MOTIF_all_c_alpha_results[compare_key]],[differ_soft,differ_all_c_alphat]
step_2_way_1_one_residue_groups,differ_group_1= compare_key_output(ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,'step_2_way_1_one_residue_group')
step_2_way_1_four_distance_residue_groups,differ_group_2= compare_key_output(ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,'step_2_way_1_four_distance_residue_group')

#%%
step_2_way_1_three_residue_groups,differ_group_3= compare_key_output(ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,'step_2_way_1_three_residue_group')
step_2_way_1_four_residue_groups,differ_group_4= compare_key_output(ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,'step_2_way_1_four_residue_group')
step_2_way_1_five_residue_groups,differ_group_5= compare_key_output(ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,'step_2_way_1_five_residue_group')
#%%
'''
Load the statistics and finalise the results
summerize for report 
'''
import pickle
import os
import numpy as np

def load_step_2_way_1_results_SUMING_UP(loading_dir_MOTIF_pick):
    load_dir = ''.join([loading_dir_MOTIF_pick,"/step_2_stat_cent_way_1_doc"])

    os.chdir('/')
    os.chdir(load_dir)
    results={}   
    results["step_2_way_1_three_residue_group"]=np.sum(pickle.load(open("step_2_way_1_three_residue_group.p", "rb")),axis=1)  
    results["step_2_way_1_four_residue_group"]=np.sum(pickle.load(open("step_2_way_1_four_residue_group.p", "rb")),axis=1)    
    results["step_2_way_1_five_residue_group"]=np.sum(pickle.load(open("step_2_way_1_five_residue_group.p", "rb")),axis=1)    
        
    results["step_2_way_1_three_distance_residue_group"]=np.sum(pickle.load(open("step_2_way_1_three_distance_residue_group.p", "rb")),axis=1)    
    results["step_2_way_1_four_distance_residue_group"]=np.sum(pickle.load(open("step_2_way_1_four_distance_residue_group.p", "rb")),axis=1)    
    results["step_2_way_1_five_distance_residue_group"]=np.sum(pickle.load(open("step_2_way_1_five_distance_residue_group.p", "rb")),axis=1) 
    for i in range(0,len(results["step_2_way_1_three_distance_residue_group"])):      
        results["step_2_way_1_three_distance_residue_group"][i]=results["step_2_way_1_three_distance_residue_group"][i]/results["step_2_way_1_three_residue_group"][i]
        results["step_2_way_1_four_distance_residue_group"][i]=results["step_2_way_1_four_distance_residue_group"][i]/results["step_2_way_1_four_residue_group"][i]
        results["step_2_way_1_five_distance_residue_group"][i]=results["step_2_way_1_five_distance_residue_group"][i]/results["step_2_way_1_five_residue_group"][i]    
    return results

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/ONGO_MOTIF_pickles_V91'
ONGO_MOTIF_soft_results=load_step_2_way_1_results_SUMING_UP(loading_dir_MOTIF_pick)

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/Surface_soft_thresh/TSG_MOTIF_pickles_V91'
TSG_MOTIF_soft_results=load_step_2_way_1_results_SUMING_UP(loading_dir_MOTIF_pick)

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/ONGO_MOTIF_pickles_V91'
ONGO_MOTIF_all_c_alpha_results=load_step_2_way_1_results_SUMING_UP(loading_dir_MOTIF_pick)

loading_dir_MOTIF_pick='C:/Users/nishy/Documents/Projects_UL/Mot_if_study_python/Summer_2020/results_and_pickles/All_C_alpha_MOTIF/TSG_MOTIF_pickles_V91'
TSG_MOTIF_all_c_alpha_results=load_step_2_way_1_results_SUMING_UP(loading_dir_MOTIF_pick)

#%% summerise for table represenattion

def update_results(key,ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,results):
    results[key][:,0]=ONGO_MOTIF_soft_results[key]
    results[key][:,1]=TSG_MOTIF_soft_results[key]
    results[key][:,2]=ONGO_MOTIF_all_c_alpha_results[key]
    results[key][:,3]=TSG_MOTIF_all_c_alpha_results[key]
    return results

results={}
results = update_results("step_2_way_1_three_residue_group",ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,results)
results = update_results("step_2_way_1_four_residue_group",ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,results)
results = update_results("step_2_way_1_five_residue_group",ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,results)
results = update_results("step_2_way_1_three_distance_residue_group",ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,results)
results = update_results("step_2_way_1_four_distance_residue_group",ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,results)
results = update_results("step_2_way_1_five_distance_residue_group",ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,results)

#%%

def compare_key_output_prob_ONGO(ONGO_MOTIF_results,TSG_MOTIF_results):
    
    differ_soft=ONGO_MOTIF_soft_results[compare_key]/(ONGO_MOTIF_soft_results[compare_key]+TSG_MOTIF_soft_results[compare_key])
    differ_all_c_alphat=ONGO_MOTIF_all_c_alpha_results[compare_key]/(ONGO_MOTIF_all_c_alpha_results[compare_key]+TSG_MOTIF_all_c_alpha_results[compare_key])

    return differ_soft,differ_all_c_alphat
prob_differ_group_3= compare_key_output_prob_ONGO(ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,'step_2_way_1_three_residue_group')
prob_differ_group_4= compare_key_output_prob_ONGO(ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,'step_2_way_1_four_residue_group')
prob_differ_group_5= compare_key_output_prob_ONGO(ONGO_MOTIF_soft_results,TSG_MOTIF_soft_results,ONGO_MOTIF_all_c_alpha_results,TSG_MOTIF_all_c_alpha_results,'step_2_way_1_five_residue_group')
#%%
def constrain_selector_way_1(ONGO_MOTIF_results_t,TSG_MOTIF_results_t):
    '''
    This function will go through the results and find the constraint
    '''
    ongo_prob=np.zeros((21,21))
    for i in range(0,21):
        for j in range(0,21):
            if ONGO_MOTIF_results_t[i][j]+TSG_MOTIF_results_t[i][j]!=0:
                ongo_prob[i][j]=ONGO_MOTIF_results_t[i][j]/(ONGO_MOTIF_results_t[i][j]+TSG_MOTIF_results_t[i][j])
            else:
                ongo_prob[i][j]=-1
                
            if (ONGO_MOTIF_results_t[i][j]/(np.sum(ONGO_MOTIF_results_t)+np.sum(TSG_MOTIF_results_t)))>=0.1:
                if ongo_prob[i][j]>0.6 or 0<ongo_prob[i][j]<0.4:
                    print("group found with center ", residues_for_index[i], ' and occurance ' , residues_for_index[j])
                    print("ONGO probability of that kind: ",ongo_prob[i][j])
    
print("Chekcing group with 3 soft thresh")
constrain_selector_way_1(step_2_way_1_three_residue_groups[0],step_2_way_1_three_residue_groups[1])
print("Chekcing group with 3 All")
constrain_selector_way_1(step_2_way_1_three_residue_groups[2],step_2_way_1_three_residue_groups[3])

print("Chekcing group with 4 soft thresh")
constrain_selector_way_1(step_2_way_1_four_residue_groups[0],step_2_way_1_four_residue_groups[1])
print("Chekcing group with 4 All")
constrain_selector_way_1(step_2_way_1_four_residue_groups[2],step_2_way_1_four_residue_groups[3])

print("Chekcing group with 5 soft thresh")
constrain_selector_way_1(step_2_way_1_five_residue_groups[0],step_2_way_1_five_residue_groups[1])
print("Chekcing group with 5 All")
constrain_selector_way_1(step_2_way_1_five_residue_groups[2],step_2_way_1_five_residue_groups[3])