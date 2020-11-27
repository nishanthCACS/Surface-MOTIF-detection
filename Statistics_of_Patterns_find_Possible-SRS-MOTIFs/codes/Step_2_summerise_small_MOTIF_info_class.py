# -*- coding: utf-8 -*-
"""
Created on %%(28-June-2020) at 03.34 P.m

@author: %A.Nishanth C00294860
"""

import os
import pickle
from copy import deepcopy
import numpy as np

class Step_2_summerise_MOTIF_info:
    """
    USE MSMS tool(1996) to calculte the soluable access area of the surface
        depth of the C_alpha carbons, and 
              of the residue 
                          from the surface
    this vertion fix the principle direction
    """
    def __init__(self, loading_dir_MOTIF_pick, saving_dir, pdb_name,consider_length_MOTIF=5):
        self.saving_dir = saving_dir
        self.pdb_name  = pdb_name #name of the class
       
        consider_length_MOTIF=consider_length_MOTIF
        '''Assign property way'''
        os.chdir('/')
        os.chdir(loading_dir_MOTIF_pick)
        MOTIF_tilted_coordinates =pickle.load(open(''.join(["MOTIF_tilted_coordinates_",self.pdb_name,".p"]), "rb"))  
        MOTIF_res =pickle.load(open(''.join(["MOTIF_res_",self.pdb_name,".p"]), "rb"))  
        MOTIF_quarter =pickle.load(open(''.join(["MOTIF_quarter_",self.pdb_name,".p"]), "rb"))  
        MOTIF_group_info =pickle.load(open(''.join(["MOTIF_group_info_",self.pdb_name,".p"]), "rb"))  
        #% take one of the MOTIF group
        group_count_info=np.zeros((21)) 
 
        PDB_MOTIF_overall_info_step_2={}     
        MOTIF_groups_avail=list(set(MOTIF_group_info))
        left_groups=0
        left_number_residues=0
        total_MOTIF_groups_residues_considered=0
        total_MOTIF_groups_considered=0
        MOTIF_step_2_cosidered=[]
        for k in range(0,len(MOTIF_groups_avail)):
            taken_MOTIF_group=MOTIF_groups_avail[k]  
            count=0
            for m in MOTIF_group_info:
                if m==taken_MOTIF_group:
                    count=count+1
            if count<21: 
                #the group should have atleat one residue
#                group_count_info[count-1]=group_count_info[count]+1
                group_count_info[count-1]=group_count_info[count-1]+1#fixing bug on 16-Oct-2020

            else:
                group_count_info[20]= group_count_info[20]+1
                
            if count<=consider_length_MOTIF:
                j=0
                taken_MOTIF_summery={}
                taken_MOTIF_group_coordinates=np.empty((count,3))
                taken_MOTIF_group_aminoacids=[]
                taken_MOTIF_group_quater_info=np.empty((count))
                for i in range(0,len(MOTIF_group_info)):
                    if MOTIF_group_info[i]==taken_MOTIF_group:
                        taken_MOTIF_group_coordinates[j][:]=MOTIF_tilted_coordinates[i,:]
                        taken_MOTIF_group_aminoacids.append(MOTIF_res[i])
                        taken_MOTIF_group_quater_info[j]= self.helper_find_quater(MOTIF_quarter,i)
                        j=j+1
                total_MOTIF_groups_residues_considered=total_MOTIF_groups_residues_considered+count
                total_MOTIF_groups_considered=total_MOTIF_groups_considered+1
                
                #% find the center of gravity of the MOTIF group and its corresponding aminoacid     
                '''
                consider/avoid the MOIF-group depend on their distance of the quarter
                First take the center of gravity and find the nearest residue
                
                check the group fell into the 
                '''
                cent_gravity=np.sum(taken_MOTIF_group_coordinates,axis=0)/8
                dist=np.empty((count))
                for j in range(0,len(taken_MOTIF_group_coordinates)):
                    dist[j]=np.linalg.norm(taken_MOTIF_group_coordinates[j]-cent_gravity)
                #    print(np.linalg.norm(taken_MOTIF_group_coordinates[j]-cent_gravity))
                #    print(np.sqrt((taken_MOTIF_group_coordinates[j][0]-cent_gravity[0])**2+ (taken_MOTIF_group_coordinates[j][1]-cent_gravity[1])**2 + (taken_MOTIF_group_coordinates[j][2]-cent_gravity[2])**2))
                
                #% define center of gravity based on the number of residues in the MOTIF group
                center_residue = np.argmin(dist)
                dist_center_residue=np.empty((count))
                for j in range(0,len(taken_MOTIF_group_coordinates)):    
                    dist_center_residue[j]=np.linalg.norm(taken_MOTIF_group_coordinates[j]-taken_MOTIF_group_coordinates[center_residue])
            
                taken_MOTIF_summery['dist_center_residue']=dist_center_residue
                taken_MOTIF_summery['taken_MOTIF_group_quater_info']=taken_MOTIF_group_quater_info
                taken_MOTIF_summery['taken_MOTIF_group_aminoacids']=taken_MOTIF_group_aminoacids
                MOTIF_step_2_cosidered.append(deepcopy(taken_MOTIF_summery))
            else:
                left_groups=left_groups+1
                left_number_residues=left_number_residues+count
                
        # the MOTIF groups center of gravity calculated
        PDB_MOTIF_overall_info_step_2['MOTIF_step_2_cosidered']=MOTIF_step_2_cosidered
        # the MOTIF groups considered for calculation
        PDB_MOTIF_overall_info_step_2['total_MOTIF_groups_considered']=total_MOTIF_groups_considered
        # the MOTIF groups considered for calculation's number residues there
        PDB_MOTIF_overall_info_step_2['total_MOTIF_groups_residues_considered']=total_MOTIF_groups_residues_considered
        # the MOTIF groups left in calculation due to they have more residues than consider_length_MOTIF
        PDB_MOTIF_overall_info_step_2['left_groups']=left_groups
        PDB_MOTIF_overall_info_step_2['left_number_residues']=left_number_residues
        # the information of overall MOTIF structures and their residues count according to the index

        # group_count_info[0]= # of the 1 residue MOTIF groups
        # group_count_info[1]= # of the 2 residue MOTIF groups
        #     :
        # group_count_info[19]= # of the 20 residue MOTIF groups
        # group_count_info[20]= # of the residues MOTIF groups > 20
        PDB_MOTIF_overall_info_step_2['group_count_info']=group_count_info

        os.chdir('/')
        if not os.path.isdir(saving_dir):
            os.makedirs(saving_dir)
        os.chdir('/')
        os.chdir(saving_dir)

        pickle.dump(group_count_info, open( ''.join(["group_count_info_",self.pdb_name,".p"]), "wb" )) 
        pickle.dump(PDB_MOTIF_overall_info_step_2, open( ''.join(["PDB_MOTIF_overall_info_step_2_",self.pdb_name,".p"]), "wb" )) 

    
    def helper_find_quater(self,MOTIF_quarter,i):
        for q in range(0,8):
           if MOTIF_quarter[q,i]==1:
               return 
           
class Step_2_group_samples_MOTIF_info:
    """
    USE MSMS tool(1996) to calculte the soluable access area of the surface
        depth of the C_alpha carbons, and 
              of the residue 
                          from the surface
    this vertion fix the principle direction
    """
    def __init__(self, loading_dir_MOTIF_pick, saving_dir, pdb_name,consider_length_MOTIF=4):
        self.saving_dir = saving_dir
        self.pdb_name  = pdb_name #name of the class            
        self.break_cons=False
     
        consider_length_MOTIF=consider_length_MOTIF
        '''Assign property way'''
        os.chdir('/')
        os.chdir(loading_dir_MOTIF_pick)
        MOTIF_group_info =pickle.load(open(''.join(["MOTIF_group_info_",self.pdb_name,".p"]), "rb"))  
        #% take one of the MOTIF group
        PDB_MOTIF_overall_info_step_2={}     
        MOTIF_groups_avail=list(set(sum(MOTIF_group_info,[])))

        for k in range(0,len(MOTIF_groups_avail)):
            taken_MOTIF_group=MOTIF_groups_avail[k]  
            count=0
            for chk_group_info in MOTIF_group_info:
                for m in chk_group_info:
                    if m==taken_MOTIF_group:
                        count=count+1
            if count==consider_length_MOTIF:
                os.chdir('/')
                os.chdir(saving_dir)
                PDB_MOTIF_overall_info_step_2 = pickle.load(open( ''.join(["PDB_MOTIF_overall_info_step_2_",self.pdb_name,".p"]), "rb" )) 
                self.PDB_MOTIF_overall_info_step_2=PDB_MOTIF_overall_info_step_2
                for taken_MOTIF_summery in PDB_MOTIF_overall_info_step_2['MOTIF_step_2_cosidered']:
                    dist_cent_res= taken_MOTIF_summery['dist_center_residue']
                    for i in range(0,len(dist_cent_res)):
                        if dist_cent_res[i]==0.0:
                            if taken_MOTIF_summery['taken_MOTIF_group_aminoacids'][i]=="C" and len(taken_MOTIF_summery['taken_MOTIF_group_aminoacids'])==consider_length_MOTIF:
                                if "Q" in  taken_MOTIF_summery['taken_MOTIF_group_aminoacids']:
                                    self.break_cons=True
                                    print(taken_MOTIF_summery['taken_MOTIF_group_aminoacids'])
        if self.break_cons:
            print(pdb_name," found")
