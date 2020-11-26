# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %A.Nishanth C00294860
"""

from motif_group_intialize_pikle import motif_group_intialize_pikle
import os
import pickle
import numpy as np
import copy#to make deep copy
#%% to run the script
working_dir = "C:/Users/nishy/Documents/Projects_UL/mot_if_STUDY/TSG"
saving_dir = "C:/Users/nishy/Documents/Projects_UL/mot_if_STUDY/TSG/finalised_pikle_results_groups" 

os.chdir('/')
os.chdir(working_dir)
print(os.getcwd())
# first check the files in the directory
coordinate_files=[]
for l in os.listdir():
    if l.endswith("_surf_atoms.csv"):
        coordinate_files.append(l)
'''
To selecting the file names
'''
for a in range(0,len(coordinate_files)):
#a=0
    os.chdir('/')
    os.chdir(working_dir)
    name_cor=coordinate_files[a]#to open the coordinate file
    #%check whether all names are same pdb id
    id_name_1=name_cor.split("_")
    pdb_name = id_name_1[0]
    print("---------------------///////--------progress number:  ",a)
    print("")
    print("PDB id in progress: ",id_name_1[0])
    # creating the object
    pdb_obj = motif_group_intialize_pikle(working_dir,saving_dir,pdb_name)
    #create the surface by triangle
    pdb_obj.creating_surface_by_triangle()
    pdb_obj.group_amino_acid_pikle()
##%% load the pijkle for see the results
##os.chdir('/')
##os.chdir(saving_dir)
##name1 = ''.join([pdb_name,"_group_",str(1),".p"])
##laoded_group = pickle.load(open(name1, "rb"))
#""""
#name convention:
#    group -is the set elements made by neighbours
#    set - represents the PDB groups
#
# then take the list from the sets and cluster them according to their length
# so the minimum number in the group allowed is 3 , 
# take the group first and go through all pikles(sets) and find find howmany time hit in each pdbs 
# -with that keep a checked list which cut that group checking again
# -When the search is going in next PDB it only check the following pikles only
#
#> maintain a check list with set of groups in the PDBs
#
#""" 
##%
#def checking_the_groups_same(group_1,group_2):
#    """
#    This fuiction take two groups and chekc whether those groups are same or not
#    If the groups are smae return 1
#    else return 0
#    
#    First take the group and the means of the groups are same
#    If the the mean is same only it goes through the elements and find out are they same
#    
#    """
#    checking = []
#    if np.mean(group_1) == np.mean(group_2):
#        if len(group_1) == len(group_2):
#            for i in range(0,len(group_1)):
#                if group_1[i]==group_2[i]:
#                    checking.append(1)
#                else:
#                    return 0
#        else:
#            return 0
#    else:
#        return 0
#    if len(checking) == len(group_1):
#        return 1
#
#
#
##%
#highest_hit = 0
#for a in range(4,20):
#    #a=0
#    pikle_files=[]
#    for l in os.listdir():
#        if l.endswith(''.join(["_group_",str(a),".p"])):
#            pikle_files.append(l)
#    PDB_set_info = [] # this has the information about the PDB sets checked 
#                        # each list represent the different PDB ids (different set) and the corresponding index 
#                        #in the array gives the information about that set
#                        #is already satisfied not unique(earlier checked with other PDBs)       
#   
#    # make the arrays for keeping the information checked group information
#    # this check the set including itself with other remaining PDB sets
#    for j in range(0,len(pikle_files)):
#        PDB_set_taken = pickle.load(open(pikle_files[j], "rb"))
#        temp_primary_set_chk = np.zeros((1,len(PDB_set_taken)))
#        PDB_set_info.extend(copy.deepcopy(temp_primary_set_chk))
#    m=0
#    p=0
#    for i in range(0,len(pikle_files)):
#        #this is containing the PDB set started for checking(primary groups as unique)
#        primary_PDB_set_checking = pickle.load(open(pikle_files[i], "rb"))
#    
#        # then go through the check list of the primary group list(primary list = set) and
#        # find out which are unique and go and check them
#        for k in range(0,len(PDB_set_info[i])):
#            if PDB_set_info[i][k] == 0:                           
#                # first load the group this is the unique group gonna check
#                primary_loaded_group = primary_PDB_set_checking[k][1]    
#                Total_count_group = 0 #to store the total occcurance of the group
#                #beacuse if primary_chk is 1 then the group is already checked
#                staisfied_PDB_ids = []
#                for j in range(i,len(pikle_files)):
#                    PDB_set_taken = pickle.load(open(pikle_files[j], "rb"))
#                    count_group = 0 # to count the how many occurance in the checking PDB set
#                    # check whether the group is checked or not if not then load it
#                    for n in range(0,len(PDB_set_taken)):
#                        if PDB_set_info[j][n] == 0:
#                            #load the set(groups) from the PDB file gonna check
#                            group_2 = PDB_set_taken[n][1]
#                            chk = checking_the_groups_same(primary_loaded_group,group_2)
#                            if chk == 1:
#                                PDB_set_info[j][n] = 1
#                                count_group = count_group + 1
#                                staisfied_PDB_ids.extend(pikle_files[j])
#                    Total_count_group = Total_count_group + count_group 
#                
#    
#
#    
#                if Total_count_group > 15:
#                    checked_group = copy.deepcopy(primary_loaded_group)# this contain the checked group and their corresponding groups have checked 
#                    # this containing the overall hit poits and the hits with which PDB_s and group information
#                    # this done for only contain one unique this unique num ber is started by m 
#                    pickle.dump(PDB_set_info , open(''.join(["amino_acid_",str(a),"PDB_set_info.p"]), "wb")) 
#                    p = p+1
#                    checked_group.append(Total_count_group)
#                    checked_group.append(p)
#                    checked_group.append(staisfied_PDB_ids)
#                    #the over all hits are represented in the pikle file name
#                    m = m + 1
#                    pickle.dump(checked_group , open(''.join(["amino_acid_",str(a),"_group_",str(m),"_",str(Total_count_group),"_total_hits",".p"]), "wb")) 
#                    print("progress unique morethan: ",m)
#                    if highest_hit < Total_count_group:
#                        highest_hit = Total_count_group
#                    print("Amino acid checking: ", a, " Coressponding i: ",str(i), " highest_hit achived: ", str(highest_hit))
#                
##            os.chdir(self.saving_dir)    
#
#
#
##%%
#"""
#This under script not working due to the PDB_set_info
#
#"""
##PDB_set_info = [] # this has the information about the PDB sets checked 
##                    # each list represent the different PDB ids (different set) and the corresponding index 
##                    #in the array gives the information about that set
##                    #is already satisfied not unique(earlier checked with other PDBs)          
##m=0
##
##for i in range(0,len(pikle_files)):
##    #this is containing the PDB set started for checking(primary groups as unique)
##    primary_PDB_set_checking = pickle.load(open(pikle_files[i], "rb"))
##
##    # make the arrays for keeping the information checked group information
##    if i == 1:
##        # this check the set including itself with other remaining PDB sets
##        for j in range(0,len(pikle_files)):
##            PDB_set_taken = pickle.load(open(pikle_files[j], "rb"))
##            temp_primary_set_chk = np.zeros((1,len(PDB_set_taken)))
##            PDB_set_info.extend(copy.deepcopy(temp_primary_set_chk))
##
##    # then go through the check list of the primary group list(primary list = set) and
##    # find out which are unique and go and check them
##    for k in range(0,len(PDB_set_info[i])):
##        if PDB_set_info[i][k] == 0:                           
##            # first load the group this is the unique group gonna check
##            primary_loaded_group = primary_PDB_set_checking[k][1]    
##            Total_count_group = 0 #to store the total occcurance of the group
##            #beacuse if primary_chk is 1 then the group is already checked
##            for j in range(i,len(pikle_files)):
##                PDB_set_taken = pickle.load(open(pikle_files[j], "rb"))
##                count_group = 0 # to count the how many occurance in the checking PDB set
##                # check whether the group is checked or not if not then load it
##                for n in range(0,len(PDB_set_taken)):
##                    if PDB_set_info[j][n] == 0:
##                        #load the set(groups) from the PDB file gonna check
##                        group_2 = PDB_set_taken[n][1]
##                        chk = checking_the_groups_same(primary_loaded_group,group_2)
##                        if chk == 1:
##                            PDB_set_info[j][n] = 1
##                            count_group = count_group + 1 
##                        PDB_set_info[j][n] = 1
##                        count_group = count_group + 1 
##                Total_count_group = Total_count_group + count_group 
##            
##            checked_group = copy.deepcopy(primary_loaded_group)# this contain the checked group and their corresponding groups have checked 
##            # this containing the overall hit poits and the hits with which PDB_s and group information
##            # this done for only contain one unique this unique num ber is started by m 
##            #the over all hits are represented in the pikle file name
##            m = m + 1
##            checked_group.append(Total_count_group)
##            checked_group.append(m)
##            os.chdir('/')
##            os.chdir(saving_dir)
##            pickle.dump(checked_group , open(''.join(["amino_acid_",str(a),"_group_",str(m),"_",str(Total_count_group),"_total_hits",".p"]), "wb")) 
###            os.chdir(self.saving_dir)