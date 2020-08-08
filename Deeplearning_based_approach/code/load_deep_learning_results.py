# -*- coding: utf-8 -*-
"""
Created on %25 Apr 2018(9:41pm)

@author: %A.Nishanth C00294860
"""
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

working_dir = "saving directory.../Deep_learning_data_initialization"
os.chdir('/')
os.chdir(working_dir)

#%% load the results and check the performance
final_batch_index = pickle.load(open("final_batch_index.p", "rb")) 
final_batch_corresponding_prediction = pickle.load(open("final_batch_corresponding_prediction.p", "rb"))  
accuracy_vector_train = pickle.load(open("accuracy_vector_train.p", "rb"))  
loss_vector = pickle.load(open("loss_vector.p", "rb"))  
#%%
plt.plot(accuracy_vector_train)
plt.ylabel("Accuarcy")
plt.xlabel("50_iteration_as_steps")
plt.show()
#%%
plt.plot(loss_vector )
plt.ylabel("Loss")
plt.xlabel("50_iteration_as_steps")
plt.show()
