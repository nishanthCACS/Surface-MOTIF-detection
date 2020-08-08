# -*- coding: utf-8 -*-
"""
Created on %25 Apr 2018(12:05pm)

@author: %A.Nishanth C00294860
"""
import tensorflow as tf
import os
import pickle
import numpy as np
#%% then load the training samples subgroups from the directory
"""
Here use every think as training sample inorder to identify the MOTIFs among them

We can do in pleanty of ways

    Just choose the PDB_ids subgroup as sets
    use 70 % of the PDB_ids_groups as training
    use 20 % of the PDB_ids_groups as testing
    use 10 % of the PDB_ids_groups as validation
    
1-Way:
    
    Randomize all groups with label and 
    
    Train with the training set the
    
    >>  loss is normal subgroup missclasification

2-Way:
    
    PDBgroup set is whole one set not-seperable
    
    Train with the training set 
    Then the subgroup is classified as ONGO or TSG MOTIF  
        then the over all class(ONGO or TSG) is choosen from the mjority vote of the subgroups
        
    >>  loss is from the PDB_class speicifiaction error
    
>>> ---------------------------------------------------------------------------

The 2-way has some issues
    Since Each PDB has atleat 100 subclasses
    Training is expensive for computation
"""
saving_dir = " saving directory information .../Deep_learning_data_initialization"
"""

                    1-Way of creation is mentioned below
    
"""
def load_pikle_type(loading_dir):
    """
    Load the property file from the given directory and combine them and make a numpy array
    """
    os.chdir('/')
    os.chdir(loading_dir)
    print(os.getcwd())
    chk =  os.listdir()
    #%create the numpy array with one class
    all_type =[]
    for i in range(0,len(chk)):
        property_chk = pickle.load(open(chk[i], "rb"))
        for p in property_chk:
            all_type.append(p[0])
    return np.array(all_type)
#%%
"""
creating the data for Deep learning once it created it can be loaded from numpy array
No need to run this cell
"""
#% load the property file of ONGO
ongo_loading_dir =  "saving directory information .../onco_started/finalised_pikle_results_groups/property_results_ongo"
ongo_prop_origin = load_pikle_type(ongo_loading_dir)
#% save the numpy array of the file
os.chdir('/')
os.chdir(saving_dir)
np.save('ongo_prop_origin.npy',ongo_prop_origin)   
#% load the property file of TSG
tsg_loading_dir =  "saving directory information.../TSG/finalised_pikle_results_groups/property_results_tsg"
tsg_prop_origin = load_pikle_type(tsg_loading_dir)
#% save the numpy array
os.chdir('/')
os.chdir(saving_dir)
np.save('tsg_prop_origin.npy',tsg_prop_origin)  
#% then combine them both
all_data = np.concatenate((ongo_prop_origin, tsg_prop_origin), axis=0)
np.save('all_data.npy',all_data) 

all_label = np.zeros((len(all_data),1))
for i in range(0,len(ongo_prop_origin)):
    all_label[i][0]=1
np.save('all_label.npy',all_label) 
#%% since the label first half is ongo and the second half is TSG
saving_dir = "saving directory information.../Deep_learning_data_initialization"
os.chdir('/')
os.chdir(saving_dir)
all_data = np.load('all_data.npy') 
all_label = np.load('all_label.npy') 
#ongo_prop_origin = np.load('ongo_prop_origin.npy')
#tsg_prop_origin = np.load('tsg_prop_origin.npy')
#%% then randomize the data
permutation=np.random.permutation(all_label.shape[0])
shuffled_all_data = all_data[permutation,:]
shuffled_all_label = all_label[permutation][:]
np.save('shuffled_all_data.npy',shuffled_all_data) 
np.save('shuffled_all_label.npy',shuffled_all_label) 
np.save('permutation.npy',permutation)
