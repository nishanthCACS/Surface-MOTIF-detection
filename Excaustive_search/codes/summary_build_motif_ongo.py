# -*- coding: utf-8 -*-
"""
Created on %%01-Jan-2017

@author: %A.Nishanth C00294860
"""
import os
import pickle
import numpy as np
hashmap_edge=np.zeros((20,20))
count=0
for i in range(0,20):
  for j in range(0,20):
    if j-i-1<0:
      count = count+1
      hashmap_edge[i,j]=count
    
count=0
for j in range(0,20):
  for i in range(0,20):
      if j-i+1>0:
          count = count+1
          hashmap_edge[i,j]=count


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

def hash_to_vertex_amino(value, hashmap_edge):
    """
    This function lookup the value and findout the coressponding indices from the hashmap
    return the 2 amino acids in the ascending order of their value(this is redundant anyway inorder to look back similarity 
    I use this)
    
    hashmap_edge: contain the mapping used
    """
    for i in range(0,20):
        for j in range(0,20):
            if hashmap_edge[i][j] == value:
                if i < j:
                    return [val_amino_acid(i),val_amino_acid(j)]
                else:
                    return [val_amino_acid(j),val_amino_acid(i)]    


"""
This script is basicaly build it for checking the groups and cerating a text file with group details

since these groups are calculated with hits higher than 15

"""

#%% load the files with satisfied mot if groups
os.chdir('/')
#os.chdir('C:/Users/nishy/Documents/Projects_UL/mot_if_STUDY/onco_started/results_created_by_R')
saving_dir = "C:/Users/nishy/Documents/Projects_UL/mot_if_STUDY/onco_started/finalised_pikle_results_groups/unique_groups" 
os.chdir(saving_dir )
print(os.getcwd())

index = 0
Total_hits = 0

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
    
for a in range(0,len(pikle_files_sorted)):    
#a = 0
    os.chdir('/')
    saving_dir = "C:/Users/nishy/Documents/Projects_UL/mot_if_STUDY/onco_started/finalised_pikle_results_groups/unique_groups" 
    os.chdir(saving_dir)
    PDB_set_detail_temp = pickle.load(open(pikle_files_sorted[a], "rb"))
    # to get the which pdb ids satisfy                         
    
    """
    Given the value it have to findout the vertices and center point from that for that group
    from the grop_string name center amino acid name found
    
    The vertices are get from the indexes trace back from the hash map
    
    from the beginning of the indices to 3 indices before containing details about the group 
    """
    
    
    g_indice_initiate=[]
    for i in range(0,len(PDB_set_detail_temp)-3):#since the python index start from zero end before the length of the list
        g_indice_initiate.append(PDB_set_detail_temp[i])
        
    #% then look up from the hash map and findout the vertices
    amino_acid_group_vertices=[]
    for i in range(len(g_indice_initiate)):
        amino_acid_group_vertices.append(hash_to_vertex_amino(g_indice_initiate[i], hashmap_edge))
    
    #% last element in the PDB_set_detail_temp is the group staisfy
    
    """
    inorder to get the PDB ids satisfied the conditions
    
    """
    satisfy_group_temp=[] 
    satisfy_group_PDB_ids=[]
    group_sat_PDBs_list_temp = PDB_set_detail_temp[len(PDB_set_detail_temp)-1]
    for i in range(0,len(group_sat_PDBs_list_temp),15):
    #    print(i)
        satisfy_group_temp.append(group_sat_PDBs_list_temp[i:i+15])
        satisfy_group_PDB_ids.append(''.join(group_sat_PDBs_list_temp[i:i+4]))
        
    #%    #Create the summary list from the details build up to now

    os.chdir('/')
    os.chdir("C:/Users/nishy/Desktop/chk_python_summarry")
    
    #%
    name_list = pikle_files_sorted[a].split("_")
    hit_points =int(name_list[5])
    center_aminoacid =val_amino_acid(int(name_list[2]))
    print('create python file part ',a)
    if a==0:
        f= open(''.join(['summary_ongo_chk.txt']),"w+")
        #then create the heads for readability
        f.write('Index'+'\t'+'Hits'+'\t'+ 'Center'+ '\t'+ 'Group'+'\r')
        f.write(str(index)+'\t\t'+str(hit_points) +'\t\t' + center_aminoacid  + '\t')
        for g in amino_acid_group_vertices:
            #insert the group details
            f.write('\t'+ ''.join(g))
        f.write('\r')
    else:
        # avoid the head details
        f.write(str(index)+'\t\t'+str(hit_points) +'\t\t' + center_aminoacid  + '\t')
        for g in amino_acid_group_vertices:
            #insert the group details
            f.write('\t'+ ''.join(g))
        f.write('\r')
        
    index = index + 1 
    Total_hits = Total_hits + hit_points
f.write('\n')
f.write("Total hits: " + str(Total_hits))
f.close() 



