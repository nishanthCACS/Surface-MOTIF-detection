# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %A.Nishanth C00294860
"""
import numpy as np
import csv
import os
os.chdir('/')
os.chdir('C:/Users/nishy/Desktop/Documents/Fall_2020/codes/Level-3_RV_coefficient/results')
#%%

f = open('RV_2_ONGO_Vs_TSG_all_averaged.csv')
csv_f = csv.reader(f)

stack=[]
for row in csv_f:
    stack.append(row)
  

#%%  
for a in [0,1,2,3]:
#a=0    #Then read the center atom from the text file
    temp = stack[a]    
    coef=[]
    for i in range(2,len(temp)):
        try:
            if float(temp[i])<0:
                print(temp[0],' ',temp[1],' check group: ', i-2)
                print(float(temp[i]))
        except:
            break
