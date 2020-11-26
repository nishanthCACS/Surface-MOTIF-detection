# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %A.Nishanth C00294860
"""
import os
"""
This script is complete analysis of groups formed
"""
#%% load the groups and find out the unique MOTIF's among TSG and ONGO
# inorder to do this take a group from the ongo and check among other with the TSG
# load the summary and make th eset of groups
 #load the txt file and edit the 
saving_dir = "C:/Users/nishy/Desktop/chk_python_summarry"
os.chdir('/')
os.chdir(saving_dir)

#%% take group details from the row
def group_center_given(name):
    """
    This function take the center and group information from therows and give the details
    """
    file = open(name, 'r') 
    rows=[]
    for line in file: 
        rows.append(line) 
    center = []
    group =[]
    hits =[]
    for row in rows:
        group_t=[]
        comma = 0
        for i in range(0,len(row)):
            if row[i]==',':
                comma = comma + 1
                if comma == 1:
                    if row[i+3]==',':
                        hits.append(int(''.join([row[i+1],row[i+2]])))
                    elif row[i+4]==',':
                        hits.append(int(''.join([row[i+1],row[i+2],row[i+3]])))
                elif comma == 2:
                    if row[i+2]==',':
                        center.append(int(row[i+1]))
                    elif row[i+3]==',':
                        center.append(int(''.join([row[i+1],row[i+2]])))
                elif comma > 2:
                    if row[i+2]=='.':
                        group_t.append(int(row[i+1]))
                    elif row[i+3]=='.':
                        group_t.append(int(''.join([row[i+1],row[i+2]])))
                    elif row[i+4]=='.':
                        group_t.append(int(''.join([row[i+1],row[i+2],row[i+3]])))
        group.append(group_t)
    return group,center,hits

group_tsg,center_tsg,hits_tsg = group_center_given('summary_tsg_group_numbered_chk.txt')
group_ongo,center_ongo,hits_ongo = group_center_given('summary_ongo_group_numbered_chk.txt')
#%% compare the groups one by one
not_unique=[]
not_unique_ongo=[]
not_unique_tsg=[]
for i in range(0,len(group_ongo)):
    #first check the center atom is same if that only then check the group is same 
    # if the group are same index of the both group is attached
    for j in range(0,len(group_tsg)):
        if center_ongo[i]==center_tsg[j]:
            if group_ongo[i]==group_tsg[j]:
                not_unique.append([i,j])
                not_unique_ongo.append(i)
                not_unique_tsg.append(j)
#%% then remove the common motif's from both set
#create the summary of the groups for each
os.chdir('/')
os.chdir("C:/Users/nishy/Desktop/chk_python_summarry")

#% read the ongo summary
name ='summary_ongo_chk.txt'
file = open(name, 'r') 
rows=[]
for line in file: 
    rows.append(line) 
not_unique_ongo.sort()
new_index =-1
a=0
Total_hits = sum(hits_ongo)
f= open(''.join(['summary_ongo_unique.txt']),"w+")
for i in range(0,len(rows)):
    if i-1 == not_unique_ongo[a]:
        if a<len(not_unique_ongo)-1:
            Total_hits  = Total_hits - hits_ongo[a]
            a=a+1
    else:
        if i < len(rows)-1:
            if i==0:
                 #then create the heads for readability
                f.write('New'+'\t\t'+'Old'+'\r')
                f.write('Index'+'\t'+'Index'+'\t'+'Hits'+'\t'+ 'Center'+ '\t'+ 'Group'+'\r')
            else:
                f.write(str(new_index)+ '\t\t')
                f.write(rows[i])

            new_index = new_index +1
f.write("Total hits: " + str(Total_hits))
f.close() 
#%%
#% read the ongo summary
name ='summary_tsg_chk.txt'
file = open(name, 'r') 
rows=[]
for line in file: 
    rows.append(line) 
not_unique_tsg.sort()
new_index =-1
a=0
Total_hits = sum(hits_tsg)
f= open(''.join(['summary_tsg_unique.txt']),"w+")
for i in range(0,len(rows)):
    if i-1 == not_unique_tsg[a]:
        if a<len(not_unique_tsg)-1:
            Total_hits  = Total_hits - hits_tsg[a]
            a=a+1
    else:
        if i < len(rows)-1:
            if i==0:
                 #then create the heads for readability
                f.write('New'+'\t\t'+'Old'+'\r')
                f.write('Index'+'\t'+'Index'+'\t'+'Hits'+'\t'+ 'Center'+ '\t'+ 'Group'+'\r')
            else:
                f.write(str(new_index)+ '\t\t')
                f.write(rows[i])
            new_index = new_index +1
f.write("Total hits: " + str(Total_hits))
f.close() 
#%% check the unique pdb details