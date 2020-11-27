# -*- coding: utf-8 -*-
"""
Created on %%(29-June-2020) at 4.47 P.m

@author: %A.Nishanth C00294860
"""

import os
import pickle
from copy import deepcopy
#import numpy as np

class Step_2_summ_sub_MOTIF_stat_cen_res_dist_doc:
    def __init__(self, Hash_residue_to_index,loading_dir_MOTIF_pick,pdb_name,way_1=True):

        self.pdb_name  = pdb_name #name of the class
        self.Hash_residue_to_index=Hash_residue_to_index

        loading_dir =''.join([loading_dir_MOTIF_pick,"/step_2_summerize"])
        os.chdir('/')
        os.chdir(loading_dir)
        PDB_MOTIF_overall_info_step_2 =pickle.load(open(''.join(["PDB_MOTIF_overall_info_step_2_",self.pdb_name,".p"]), "rb"))  

        MOTIF_step_2_cosidered= PDB_MOTIF_overall_info_step_2['MOTIF_step_2_cosidered']
        #to avoid unkonwn residues come into play

        key_error=False
        if way_1:
            load_stat_dir = ''.join([loading_dir_MOTIF_pick,"/step_2_stat_cent_way_1_doc"])
            os.chdir('/')
            os.chdir(load_stat_dir)
            #to avoid unessery loading of the MOT group stat details      
            MOTIF_1=True
            MOTIF_2=True
            MOTIF_3=True
            MOTIF_4=True
            MOTIF_5=True
            
            for taken_MOTIF_summery in MOTIF_step_2_cosidered:
                dist_center_residue=  taken_MOTIF_summery['dist_center_residue']
                taken_MOTIF_group_aminoacids = taken_MOTIF_summery['taken_MOTIF_group_aminoacids']       
                for key in taken_MOTIF_group_aminoacids:
                    if key not in list(self.Hash_residue_to_index.keys()):
                        print(loading_dir_MOTIF_pick)
                        print("Residue ", key ," not reported thus PDB ",pdb_name," skipped")
                        key_error=True
                if not key_error:
                    if len(dist_center_residue)==1:
                        if MOTIF_1:
                            step_2_way_1_one_residue_group =pickle.load(open('step_2_way_1_one_residue_group.p', "rb"))  
                        index_1=self.Hash_residue_to_index[taken_MOTIF_group_aminoacids[0]]
                        step_2_way_1_one_residue_group[index_1]=step_2_way_1_one_residue_group[index_1]+1
                        MOTIF_1=False
                    elif len(dist_center_residue)==2:
                        # just update the both hits
                        if MOTIF_2:
                            step_2_way_1_two_residue_group =pickle.load(open('step_2_way_1_two_residue_group.p', "rb"))  
                            step_2_way_1_two_distance_residue_group =pickle.load(open('step_2_way_1_two_distance_residue_group.p', "rb"))  
                        MOTIF_2=False
            
                        idex_1=Hash_residue_to_index[taken_MOTIF_group_aminoacids[0]]
                        idex_2=Hash_residue_to_index[taken_MOTIF_group_aminoacids[1]]
                        if idex_1<idex_2:
                            step_2_way_1_two_residue_group[idex_1][idex_2]=step_2_way_1_two_residue_group[idex_1][idex_2]+1
                            if dist_center_residue[0]>0:
                                step_2_way_1_two_distance_residue_group[idex_1][idex_2]=dist_center_residue[0]
                            else:
                                step_2_way_1_two_distance_residue_group[idex_1][idex_2]=dist_center_residue[1]
                        else:
                            step_2_way_1_two_residue_group[idex_2][idex_1]=step_2_way_1_two_residue_group[idex_2][idex_1]+1
                            if dist_center_residue[0]>0:
                                step_2_way_1_two_distance_residue_group[idex_1][idex_2]=dist_center_residue[0]
                            else:
                                step_2_way_1_two_distance_residue_group[idex_1][idex_2]=dist_center_residue[1]
                    elif len(dist_center_residue)==3:
                        if MOTIF_3:
                            step_2_way_1_three_residue_group =pickle.load(open('step_2_way_1_three_residue_group.p', "rb"))  
                            step_2_way_1_three_distance_residue_group =pickle.load(open('step_2_way_1_three_distance_residue_group.p', "rb"))  
                        MOTIF_3=False
                        step_2_way_1_three_residue_group,step_2_way_1_three_distance_residue_group= self.helper_way_1_stat(taken_MOTIF_group_aminoacids,dist_center_residue,step_2_way_1_three_residue_group,step_2_way_1_three_distance_residue_group)
                    elif len(dist_center_residue)==4:
                        if MOTIF_4:
                            step_2_way_1_four_residue_group =pickle.load(open('step_2_way_1_four_residue_group.p', "rb"))  
                            step_2_way_1_four_distance_residue_group =pickle.load(open('step_2_way_1_four_distance_residue_group.p', "rb"))  
                        MOTIF_4=False
                        step_2_way_1_four_residue_group,step_2_way_1_four_distance_residue_group= self.helper_way_1_stat(taken_MOTIF_group_aminoacids,dist_center_residue,step_2_way_1_four_residue_group,step_2_way_1_four_distance_residue_group)
                    elif len(dist_center_residue)==5:
                        if MOTIF_5:
                            step_2_way_1_five_residue_group =pickle.load(open('step_2_way_1_five_residue_group.p', "rb"))  
                            step_2_way_1_five_distance_residue_group =pickle.load(open('step_2_way_1_five_distance_residue_group.p', "rb"))  
                        MOTIF_5=False
                        step_2_way_1_five_residue_group,step_2_way_1_five_distance_residue_group= self.helper_way_1_stat(taken_MOTIF_group_aminoacids,dist_center_residue,step_2_way_1_five_residue_group,step_2_way_1_five_distance_residue_group)
                    else:
                        print("some thing here", pdb_name)
            if not key_error:
                if not MOTIF_1:            
                    pickle.dump(step_2_way_1_one_residue_group, open('step_2_way_1_one_residue_group.p', "wb" )) 
                if not MOTIF_2:            
                    pickle.dump(step_2_way_1_two_residue_group, open('step_2_way_1_two_residue_group.p', "wb" )) 
                    pickle.dump(step_2_way_1_two_distance_residue_group, open('step_2_way_1_two_distance_residue_group.p', "wb" )) 
                if not MOTIF_3:            
                    pickle.dump(step_2_way_1_three_residue_group, open('step_2_way_1_three_residue_group.p', "wb" )) 
                    pickle.dump(step_2_way_1_three_distance_residue_group, open('step_2_way_1_three_distance_residue_group.p', "wb" )) 
                if not MOTIF_4:            
                    pickle.dump(step_2_way_1_four_residue_group, open('step_2_way_1_four_residue_group.p', "wb" )) 
                    pickle.dump(step_2_way_1_four_distance_residue_group, open('step_2_way_1_four_distance_residue_group.p', "wb" )) 
                if not MOTIF_5:                        
                    pickle.dump(step_2_way_1_five_residue_group, open('step_2_way_1_five_residue_group.p', "wb" ))           
                    pickle.dump(step_2_way_1_five_distance_residue_group, open('step_2_way_1_five_distance_residue_group.p', "wb" )) 

        else:
            '''
            For way-2 giving three center with two of them gether
            
            Thus first center residue as key
            each of the group as 21 x 21 array contain the hits
            make all the combination as possible
            e.g: 
                4-subgroup  contain center is "C" other residues as D, E, Q
                then
                center residue "C" is arraay is chosen 
                where DE, DQ, EQ are given hits
                
                 5-subgroup  contain center is "C" other residues as D, E, Q,H
                then
                center residue "C" is arraay is chosen 
                where DE, DQ, EQ,DH, EH, QH are given hits
            
            inorder to dothis first these center residue is poped
            then rest of the residues are sorted from lowest to highest depends on the mapping
            
            then take the lowest and group them with rest of the list ; then remove the lowest and repeat until 1 residue exist in the group
            '''
            
            load_stat_dir = ''.join([loading_dir_MOTIF_pick,"/step_2_stat_cent_way_2_doc"])
            os.chdir('/')
            os.chdir(load_stat_dir)
            for taken_MOTIF_summery in MOTIF_step_2_cosidered:   
                dist_center_residue=  taken_MOTIF_summery['dist_center_residue']
                
                if len(dist_center_residue)>2:
                    dist_center_residue=  taken_MOTIF_summery['dist_center_residue']
                    taken_MOTIF_group_aminoacids=deepcopy(taken_MOTIF_summery['taken_MOTIF_group_aminoacids'])
                    key_error=False
                    for key in taken_MOTIF_group_aminoacids:
                        if key not in list(self.Hash_residue_to_index.keys()):
                            print(loading_dir_MOTIF_pick)
                            print("Residue ", key ," not reported thus PDB ",pdb_name," skipped")
                            key_error=True
                    if not key_error:
                        center_residue= taken_MOTIF_group_aminoacids.pop(self.find_center_residue(dist_center_residue))
                        step_2_way_2_group =pickle.load(open(''.join([center_residue,'_way_2_five_residue_group.p']), "rb"))  
                
                        mapped_group=[]   
                        for res in taken_MOTIF_group_aminoacids:
                            mapped_group.append(Hash_residue_to_index[res])
                        mapped_group.sort(reverse=True)
    #                    sel_group_all=[]
                        while len(mapped_group)>1:
                            sel_res=mapped_group.pop()
    #                        sel_group_t=[]
                            for res in mapped_group:
    #                            sel_group_t.append([sel_res,res])
                                step_2_way_2_group[sel_res][res]= step_2_way_2_group[sel_res][res]+1
    #                        sel_group_all.append(sel_group_t)
                        
                        pickle.dump(step_2_way_2_group, open(''.join([center_residue,'_way_2_five_residue_group.p']), "wb" )) 

            
    def find_center_residue(self,dist_center_residue):
        for i in range(0,len(dist_center_residue)):
            if dist_center_residue[i]==0:
                return i
            
    def helper_way_1_stat(self,taken_MOTIF_group_aminoacids,dist_center_residue,step_2_way_1_any_residue_group,step_2_way_1_any_residue_dis_group):
        center_res_index=self.find_center_residue(dist_center_residue)
        index_1=self.Hash_residue_to_index[taken_MOTIF_group_aminoacids[center_res_index]]
    
        for k in range(0,len(dist_center_residue)):
            if k !=center_res_index:
                index_2=self.Hash_residue_to_index[taken_MOTIF_group_aminoacids[k]]
                step_2_way_1_any_residue_group[index_1][index_2]= step_2_way_1_any_residue_group[index_1][index_2]+1
                step_2_way_1_any_residue_dis_group[index_1][index_2]= step_2_way_1_any_residue_dis_group[index_1][index_2]+dist_center_residue[k]
        return step_2_way_1_any_residue_group,step_2_way_1_any_residue_dis_group