# -*- coding: utf-8 -*-
"""
Created on %28-Dec-2017

@author: %A.Nishanth C00294860
"""

import csv
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
import copy#to make deep copy
#%%
# to compare the two lists
class motif_group_intialize_pikle:
    """
    This class is basically written for creating the pikle of groups of
    aminoacids with neighbours deatlis with ascending order of length in the group

    """
    def __init__(self,working_dir,saving_dir,pdb_name):
        
        self.working_dir=working_dir
        self.saving_dir =saving_dir
        self.pdb_name = pdb_name   

        #create the hash map of the amino acid edges
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

        self.hashmap_edge = hashmap_edge
        os.chdir('/')
        os.chdir(working_dir)
        
    def amino_acid_val(self,amino_acid):
        """
        this function tells amino acid value for find the hash map value for edge
        """
    
        if(amino_acid=="G"):
          val=0
        elif (amino_acid=="M"):
          val=1
        elif(amino_acid=="R"):
          val=2
        elif(amino_acid=="K"):
          val=3
        elif(amino_acid=="D"):
          val=4
        elif(amino_acid=="E"):
          val=5
        elif(amino_acid=="Q"):
          val=6
        elif(amino_acid=="N"):
          val=7
        elif(amino_acid=="H"):
          val=8
        elif(amino_acid=="S"):
          val=9
        elif(amino_acid=="T"):
          val=10
        elif(amino_acid=="Y"):
          val=11
        elif(amino_acid=="C"):
          val=12
        elif(amino_acid=="W"):
          val=13
        elif(amino_acid=="A"):
          val=14
        elif(amino_acid=="I"):
          val=15
        elif(amino_acid=="L"):
          val=16
        elif(amino_acid=="F"):
          val=17
        elif(amino_acid=="V"):
          val=18
        elif(amino_acid=="P"):
          val=19
        return(val)  
        
        #% change the directory to load the coordinate file and the corresponding aminoacids


    #% load the coordinate file 
    def coordinate_numpy(self,name):
        """
        This function take file with name and return numpy array of that
        """
        with open(name, newline='') as f:
            rows=[]
            reader = csv.reader(f)
            for row in reader:
                rows.append(row)
        #% then store the coordinate file as np array
        coordinates_t = np.zeros((len(rows)-1,3))
        for i in range(1,len(rows)):
            coordinates_t[i-1][:]=rows[i][:]
        return coordinates_t     
    
    def aminoacid_function_retrieve(self,name_amino_acid):
        """
        This function return the surface amino acids
        """
        with open(name_amino_acid, newline='') as f:
            aminoacid_surface_t=[]
            reader = csv.reader(f)
            for row in reader:
                aminoacid_surface_t.extend(row)
            aminoacid_surface = aminoacid_surface_t[1:len(aminoacid_surface_t)]
        return aminoacid_surface
    
    #% function for find the equation of triangle
    def distance_surface(self,coordinates,triangle_points,missed_point_index):
        """
        This function calculates the diatance from the missed point
        
        coordinates: cartesian coordinates of the surface atoms
        verrtices  : Triangles points
        verrtices_position_triangle: Which triangle from the vertices wanna be considered
        missed_point_index: The point which need to calculate the distance
        """
        missed_point = coordinates[missed_point_index]
        
        tri_point_1 = coordinates[triangle_points[0]]
        tri_point_2 = coordinates[triangle_points[1]]
        tri_point_3 = coordinates[triangle_points[2]]
        
        #calculating the vectors for surface
        x = tri_point_1 - tri_point_2
        y = tri_point_1 - tri_point_3
            
        eq_t=np.cross(x,y)
        distance_t = (eq_t[0]*(missed_point[0]-tri_point_1[0])\
        + eq_t[1]*(missed_point[1]-tri_point_1[1])\
        + eq_t[2]*(missed_point[2]-tri_point_1[2]))/np.linalg.norm(eq_t)
        
        if distance_t<0:
            return (-1)*distance_t
        else:
            return distance_t


    def creating_surface_by_triangle(self):
        """
        This function create the suraface with triangles from
        1> create the convex hull using the library
        2> Extend the surface traiangulation wich should cover all surface points
        """
        coordinates =  self.coordinate_numpy(''.join([self.pdb_name,"_surf_atoms.csv"]))   
        coordinates_polar = self.coordinate_numpy(''.join([self.pdb_name,"_surf_polar_atoms.csv"]))  # the surface atoms given by R they hold the first coloumn radius 
                                                          # second coloyumn xy angle(theta) and then xy_z angle(phi)
                                          
        #% creating the convex hull using the libraray
        tri =  ConvexHull(coordinates)
        
        #% define polyhyderen methods
        verrtices=[]
        chk_points=[]
        for vert in tri.simplices:
            verrtices.append(vert)
            for point in vert:
                chk_points.append(point)
        
        #% to findout the points missed by conves hull creation
        chk_points=list(set(chk_points))
        chk_points.sort()
        missing=[]
        k=0
        for i in range(0,len(coordinates)):
            if i==chk_points[k]:
                k=k+1
                if k == len(chk_points):
                    k=k-1
            else:
                missing.append(i)
        #% then create the triangles
        #from the missing points order them as according to their radius
        missing_points_radius=[]
        for i in missing:
            missing_points_radius.append(coordinates_polar[i][0])
        #    missing_points_radius.append(    coordinates_polar[0])
        assending_order_of_radius=sorted(range(len(missing_points_radius)), key=lambda k: missing_points_radius[k])
        #% from the missing point radiusus findout the mising points 
        """
        checking_points_order: keep the missing points checking order in the decending order of the radius
        """
        checking_points_order=[]
        for i in range(0,len(assending_order_of_radius)):
            checking_points_order.append(missing[assending_order_of_radius[len(assending_order_of_radius)-i-1]])
        ##%% sorting the list
        #s = [8, 9, 1, 2, 3]
        #sorted(range(len(s)), key=lambda k: s[k])
        #% place the triangles phi angles with the triangles respectively
        """
        inorder to dothis the triangles points are selected one by one and choose the miimum and maximum phi angle
        and place with the triangle position{as hashmap triangle position to phi angles(xy_z) minimum and maximum}
        
        """
        phi_t=[]
        # first select the phe angle of the three points in triangle and place in the hashmap
        for i in range(0,len(verrtices)):
            temp=[]
            temp.append(coordinates_polar[verrtices[i][0]][2])
            temp.append(coordinates_polar[verrtices[i][1]][2])
            temp.append(coordinates_polar[verrtices[i][2]][2])
            temp.sort()
            phi_t.append(temp)
        
        # sort the theta same as well
        theta_t=[]
        for i in range(0,len(verrtices)):
            temp=[]
            temp.append(coordinates_polar[verrtices[i][0]][1])
            temp.append(coordinates_polar[verrtices[i][1]][1])
            temp.append(coordinates_polar[verrtices[i][2]][1])
            temp.sort()
            theta_t.append(temp)

        #% choose the triangles to check the ditance of the missed point
        """
         inorder to find which traingle distances have to be checked 
         - first choose the triangles have theta and phi in the range of of the missing coordinate
         - if that is more than one 
         - then find the equation of the triangles 
         - calculate the distance from the point to triangle and choose the minimum distance triangle and break that
         
         for checking purpose use uncomment and use triangle details
        
                """
        # first run through the phi angles and find out the triangle satisify the first condition
        for s in range(0,len(missing)):
            #s = 0#indices of the point in the missing
            selected_point_phi_angle = coordinates_polar[checking_points_order[s]][2]
            selected_point_theta_angle =  coordinates_polar[checking_points_order[s]][1]
            triangle_sat=[]# that triangle position in the list is added
            #triangle_details=[]
            for i in range(0,len(phi_t)):
                if phi_t[i][0] <= selected_point_phi_angle <= phi_t[i][2]:
                    if theta_t[i][0] <= selected_point_theta_angle <= theta_t[i][2]:
                        triangle_sat.append(i)
            #            triangle_details.append([i,phi_t[i][0],phi_t[i][2], theta_t[i][0],theta_t[i][2]] )
#            """
#            If the satisfied triangle in these condition not staisfied loose the constraints
#            
#            choose only the phi angle and theta angle seperately
#            """
#
#            if len(triangle_sat)==0:
#                for i in range(0,len(phi_t)):
#                    if phi_t[i][0] < selected_point_phi_angle < phi_t[i][2]:
#                        triangle_sat.append(i)
#                    if theta_t[i][0] < selected_point_theta_angle < theta_t[i][2]:
#                        triangle_sat.append(i)
            """
            Even if that case not satisfied loose the dondition even further
            Here only consider the 15 degree in all ways phi +15 degree and -15 degree
                                                        theata +15 degree and -15 degree
            check which triangles middle angle fell into these ranges 
            
            25 degree = 0.2618 radians
                    """  
            if len(triangle_sat)==0:
        #        print("PDB-- 30 degree condition: ",pdb_name)
        #        print("corfessponding i:::::::::: ",i)
                for i in range(0,len(phi_t)):
                    phi_condition=0
                    for o in range(0,3):
                        if selected_point_phi_angle - 0.2618 <= phi_t[i][o] <= selected_point_phi_angle + 0.2618:
                            phi_condition=1
                    if phi_condition == 1:
                        for o in range(0,3):
                            if selected_point_theta_angle - 0.2618  <= theta_t[i][o] <= selected_point_theta_angle + 0.2618:
                                triangle_sat.append(i)
                                
        
            if len(triangle_sat)==0:
        #        print("PDB-- 60 degree condition: ",pdb_name)
        #        print("corfessponding i:::::::::: ",i)
                for i in range(0,len(phi_t)):
                    phi_condition=0
                    for o in range(0,3):
                        if selected_point_phi_angle - 0.5236  <= phi_t[i][o] <= selected_point_phi_angle + 0.5236 :
                            phi_condition=1
                    if phi_condition == 1:
                        for o in range(0,3):
                            if selected_point_theta_angle - 0.5236  <= theta_t[i][o] <= selected_point_theta_angle + 0.5236 :
                                triangle_sat.append(i)
                            
            if len(triangle_sat)==0:
        #        print("PDB-- 90 degree condition: ",pdb_name)
        #        print("corfessponding i:::::::::: ",i)
                for i in range(0,len(phi_t)):
                    phi_condition=0
                    for o in range(0,3):
                        if selected_point_phi_angle - 0.7853  <= phi_t[i][o] <= selected_point_phi_angle + 0.7853  :
                            phi_condition=1
                    if phi_condition == 1:
                        for o in range(0,3):
                            if selected_point_theta_angle - 0.7853  <= theta_t[i][o] <= selected_point_theta_angle + 0.7853  :
                                triangle_sat.append(i)
                                
            if len(triangle_sat)==0:
                print("PDB-- 120 degree condition: ",self.pdb_name)
                print("corfessponding i:::::::::: ",i)
                for i in range(0,len(phi_t)):
                    phi_condition=0
                    for o in range(0,3):
                        if selected_point_phi_angle - 1.0472 <= phi_t[i][o] <= selected_point_phi_angle + 1.0472  :
                            phi_condition=1
                    if phi_condition == 1:
                        for o in range(0,3):
                            if selected_point_theta_angle - 1.0472  <= theta_t[i][o] <= selected_point_theta_angle + 1.0472  :
                                triangle_sat.append(i) 
                                
            if len(triangle_sat)==0:
                print("PDB-- 150 degree condition: ",self.pdb_name)
                print("corfessponding i:::::::::: ",i)
                for i in range(0,len(phi_t)):
                    phi_condition=0
                    for o in range(0,3):
                        if selected_point_phi_angle - 1.309 <= phi_t[i][o] <= selected_point_phi_angle + 1.309  :
                            phi_condition=1
                    if phi_condition == 1:
                        for o in range(0,3):
                            if selected_point_theta_angle - 1.309  <= theta_t[i][o] <= selected_point_theta_angle + 1.309  :
                                triangle_sat.append(i) 
            #then use the satisfied triangles to findout the distance
            #% to findout the equations
            
            distances=[]
            for i in range(0,len(triangle_sat)):
                d_t = self.distance_surface(coordinates,verrtices[triangle_sat[i]],checking_points_order[s])
                distances.append(d_t)
                    
            #%from the distances use the minimum distance 
            """
            Consider the breaking triangle contain a,b,c as points and added point gonna be d
            then the set of new triangles are
            a, b, d
            b, c, d
            a, c, d
            
            since these triangle_brak and theta and phi angles are used again and again after updated also
            thus those are made as deep copies
            """
            assending_order_of_distance_key=sorted(range(len(distances)), key=lambda k: distances[k])
            #then choose the triangle that have mnimum distance
            #then break that triangle into three triangles
            triangle_break =  copy.deepcopy(verrtices[triangle_sat[assending_order_of_distance_key[0]]])
            # here need to make deep copy for coping the elements
            #when break the trangle parellely update the phi and theta data as well
            # using deepcopy to deep copy 
            phi_in_rem_triangle = copy.deepcopy(phi_t[triangle_sat[assending_order_of_distance_key[0]]])
            theta_in_rem_triangle =copy.deepcopy(theta_t[triangle_sat[assending_order_of_distance_key[0]]])
            #print("Phi angle in the removed place1: ",phi_in_rem_triangle)
            #%
            #placing the first triangle a,b,d
            verrtices.insert(triangle_sat[assending_order_of_distance_key[0]], np.asarray([triangle_break[0],triangle_break[1],checking_points_order[s]]))
            phi_gonna_insert = [phi_in_rem_triangle[0],phi_in_rem_triangle[1],selected_point_phi_angle]
            theta_gonna_insert = [theta_in_rem_triangle[0],theta_in_rem_triangle[1],selected_point_theta_angle]
            phi_gonna_insert.sort()
            theta_gonna_insert.sort()
            """
            checking purpose
            #print("inserted_phi_0:  ",phi_gonna_insert)
            #print("inserted_theta_0:  ",theta_gonna_insert)
            """
            phi_t.insert(triangle_sat[assending_order_of_distance_key[0]],phi_gonna_insert)
            theta_t.insert(triangle_sat[assending_order_of_distance_key[0]],theta_gonna_insert)
            """
            checking purpose
            print("inserted_phi_1:  ",phi_gonna_insert)
            print("inserted_theta_1:  ",theta_gonna_insert)
            """
            #%
            ##placing the second triangle b,c,d
            verrtices.insert(triangle_sat[assending_order_of_distance_key[0]], np.asarray([triangle_break[1],triangle_break[2],checking_points_order[s]]))
            phi_gonna_insert = [phi_in_rem_triangle[1],phi_in_rem_triangle[2],selected_point_phi_angle]
            theta_gonna_insert = [theta_in_rem_triangle[1],theta_in_rem_triangle[2],selected_point_theta_angle]
            phi_gonna_insert.sort()
            theta_gonna_insert.sort()
            phi_t.insert(triangle_sat[assending_order_of_distance_key[0]],phi_gonna_insert)
            theta_t.insert(triangle_sat[assending_order_of_distance_key[0]],theta_gonna_insert)
            """
            checking purpose
            print("inserted_phi_2:  ",phi_gonna_insert)
            print("inserted_theta_2:  ",theta_gonna_insert)
            """
            ##placing the third triangle a,c,d
            verrtices.insert(triangle_sat[assending_order_of_distance_key[0]], np.asarray([triangle_break[0],triangle_break[2],checking_points_order[s]]))
            phi_gonna_insert = [phi_in_rem_triangle[0],phi_in_rem_triangle[2],selected_point_phi_angle]
            theta_gonna_insert = [theta_in_rem_triangle[0],theta_in_rem_triangle[2],selected_point_theta_angle]
            phi_gonna_insert.sort()
            theta_gonna_insert.sort()
            phi_t.insert(triangle_sat[assending_order_of_distance_key[0]],phi_gonna_insert)
            theta_t.insert(triangle_sat[assending_order_of_distance_key[0]],theta_gonna_insert)
            """
            checking purpose
            print("inserted_phi_3:  ",phi_gonna_insert)
            print("inserted_theta_3:  ",theta_gonna_insert)
            """
            #%
            # remove the triangle and phi and theta angles
            del(verrtices[triangle_sat[assending_order_of_distance_key[0]]+3])
            del(phi_t[triangle_sat[assending_order_of_distance_key[0]]+3])
            del(theta_t[triangle_sat[assending_order_of_distance_key[0]]+3])
            """
            print("s: ",s)
            print("triangle added length: ",len(verrtices))
               
            """
        self.coordinates = coordinates      
        self.coordinates_polar = coordinates_polar      
        self.tri = tri
        self.vertices = verrtices
#        return  coordinates, coordinates_polar, tri, verrtices  


    def group_amino_acid_pikle(self):
        """
        This fuction create the pikle of the group of triangle points and center
        
        go through the surface atoms and findoutthe hit map of triangles in triangles
        """
        
        amino_acid_details = self.aminoacid_function_retrieve(''.join([self.pdb_name,"_aminoacids.csv"]))# to save the amino acid details
        
        coordinates = self.coordinates
        verrtices = self.vertices
        hashmap_edge = self.hashmap_edge
        
        cooridnate_triangle_sat_t = np.zeros((len(coordinates),len(verrtices)))
        for i in range(0,len(coordinates)):
            for j in range(0,len(verrtices)):
                # if the coordinate point in the triangle set sat = 1
                sat = 0
                for vert in verrtices[j]:
                    if vert == i:
                        sat = 1
                if sat == 1: 
                    cooridnate_triangle_sat_t[i][j] = 1 
        
        #% then from the selected triangles select the edges which doesn't have center atom as one set
        group_from_coordinates=[]
        for i in range(0,len(cooridnate_triangle_sat_t)):
            group_temp = []
            for j in range(0,len(verrtices)):
                if cooridnate_triangle_sat_t[i][j]==1:
                   temp = []
                   for positions in verrtices[j]:
                       if positions != i:
                           temp.append(positions)
                   group_temp.append(temp)
            group_from_coordinates.append(group_temp)
        
        #% then change the values to aminoacids
        group = []
        group_details=[]
        for i in range(0,len(coordinates)):
            temp=[]
            # first append the center atom amino acid details
            for chk_group in group_from_coordinates[i]:
                #take the group detail and find out the index of aminoacid
                index_1 = self.amino_acid_val(amino_acid_details[chk_group[0]])
                index_2 = self.amino_acid_val(amino_acid_details[chk_group[1]])
                temp.append(hashmap_edge[index_1][index_2])
            temp.sort()# because this only have the neighbour thus easy to findout the same group
            temp_2 = []
            temp_2.append(hashmap_edge[self.amino_acid_val(amino_acid_details[i])][self.amino_acid_val(amino_acid_details[i])])
            temp_2.append(temp)
            group.append(temp_2)
            group_details.append([self.amino_acid_val(amino_acid_details[i]),len(temp)])
        #% then cluster the groups
        #first group the center aminoacids together with ascednding order of number of elements in the group
        for j in range(0,20):    
            length_group_temp = []
            length_group_temp_index = []
            sub_group_temp = []
            for i in range(0, len(group_details)):
                if group_details[i][0]==j:
                    length_group_temp.append(group_details[i][1])
                    length_group_temp_index.append(i)
            #then sort the length of the groupa nd find the coressponding index keys and order them accordingly
            assending_order_of_group_length_key=sorted(range(len(length_group_temp)), key=lambda k: length_group_temp[k])
            # from the assending_order_of_group_length_key find the orginal group position
            for key_len in assending_order_of_group_length_key:
                sub_group_temp.append(group[length_group_temp_index[key_len]])

#            change the dirsctory to group saving directory

            os.chdir('/')
            os.chdir(self.saving_dir)
            pickle.dump(sub_group_temp, open(''.join([self.pdb_name,"_group_",str(j),".p"]), "wb")) 
            
        #% plot the resukts and save it
        def plot_surface(self):
            """
            Make this as fucntion to plot and save if needed
            plot the results and save it
            """    
            coordinates = self.coordinates
            verrtices = self.verrtices
            
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection='3d')
            
            ## The triangles in parameter space determine which x, y, z points are
            ## connected by an edge
            ax.plot_trisurf(coordinates[:,0],coordinates[:,1], coordinates[:,2], triangles=self.tri.simplices, cmap=plt.cm.Spectral)
            plt.show()
            fig.savefig(''.join([self.pdb_name,'_only_convex.png']))

#            plot after the inserted triangles convex is changed

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection='3d')
            ax.plot_trisurf(coordinates[:,0],coordinates[:,1], coordinates[:,2], triangles=verrtices, cmap=plt.cm.Spectral)
            plt.show()
            fig.savefig(''.join([self.pdb_name,'After_triangulated.png']))

class motif_group_find_pikle:
    """
    This class is basically written for check the pikle of groups of
    aminoacids PDB files
    
    And finds out the hitpoints figure the groups have highest hit prefered

    > hits_satisfy  : this basically check the howmany counts the group should mnimally have
    > amino_acid    : Which aminoacid central group gonna check that contains(0-19)
    """
    def __init__(self,working_dir,saving_dir,hits_satisfy,amino_acid_group):
        self.working_dir=working_dir
        self.saving_dir =saving_dir
        self.hits_satisfy = hits_satisfy  
        self.amino_acid_group = amino_acid_group
        print("Intialized")
        
        
    def checking_the_groups_same(self,group_1,group_2):
        """
        This fuiction take two groups and chekc whether those groups are same or not
        If the groups are smae return 1
        else return 0
        
        First take the group and the means of the groups are same
        If the the mean is same only it goes through the elements and find out are they same
        
        """
        checking = []
        if np.mean(group_1) == np.mean(group_2):
            if len(group_1) == len(group_2):
                for i in range(0,len(group_1)):
                    if group_1[i]==group_2[i]:
                        checking.append(1)
                    else:
                        return 0
            else:
                return 0
        else:
            return 0
        if len(checking) == len(group_1):
            return 1

    def find_unique_group(self):      
        
        """"
        name convention:
            group -is the set elements made by neighbours
            set - represents the PDB groups
        
         then take the list from the sets and cluster them according to their length
         so the minimum number in the group allowed is 3 , 
         take the group first and go through all pikles(sets) and find find howmany time hit in each pdbs 
         -with that keep a checked list which cut that group checking again
         -When the search is going in next PDB it only check the following pikles only
        
        > maintain a check list with set of groups in the PDBs
        
        """ 
        print("In the function")
        os.chdir('/')
        os.chdir(self.working_dir)
        highest_hit = 0
        #for a in range(4,20):
        a = self.amino_acid_group
        pikle_files=[]
        for l in os.listdir():
            if l.endswith(''.join(["_group_",str(a),".p"])):
                pikle_files.append(l)
        PDB_set_info = [] # this has the information about the PDB sets checked 
                            # each list represent the different PDB ids (different set) and the corresponding index 
                            #in the array gives the information about that set
                            #is already satisfied not unique(earlier checked with other PDBs)       
           
        # make the arrays for keeping the information checked group information
        # this check the set including itself with other remaining PDB sets
        for j in range(0,len(pikle_files)):
            PDB_set_taken = pickle.load(open(pikle_files[j], "rb"))
            temp_primary_set_chk = np.zeros((1,len(PDB_set_taken)))
            PDB_set_info.extend(copy.deepcopy(temp_primary_set_chk))
        m=0
        p=0
        for i in range(0,len(pikle_files)):
            #this is containing the PDB set started for checking(primary groups as unique)
            primary_PDB_set_checking = pickle.load(open(pikle_files[i], "rb"))
        
            # then go through the check list of the primary group list(primary list = set) and
            # find out which are unique and go and check them
            for k in range(0,len(PDB_set_info[i])):
                if PDB_set_info[i][k] == 0:                           
                    # first load the group this is the unique group gonna check
                    primary_loaded_group = primary_PDB_set_checking[k][1]    
                    Total_count_group = 0 #to store the total occcurance of the group
                    #beacuse if primary_chk is 1 then the group is already checked
                    staisfied_PDB_ids = []
                    for j in range(i,len(pikle_files)):
                        PDB_set_taken = pickle.load(open(pikle_files[j], "rb"))
                        count_group = 0 # to count the how many occurance in the checking PDB set
                        # check whether the group is checked or not if not then load it
                        for n in range(0,len(PDB_set_taken)):
                            if PDB_set_info[j][n] == 0:
                                #load the set(groups) from the PDB file gonna check
                                group_2 = PDB_set_taken[n][1]
                                chk = self.checking_the_groups_same(primary_loaded_group,group_2)
                                if chk == 1:
                                    PDB_set_info[j][n] = 1
                                    count_group = count_group + 1
                                    staisfied_PDB_ids.extend(pikle_files[j])
                        Total_count_group = Total_count_group + count_group 

                    if Total_count_group > self.hits_satisfy:
                        checked_group = copy.deepcopy(primary_loaded_group)# this contain the checked group and their corresponding groups have checked 
                        # this containing the overall hit poits and the hits with which PDB_s and group information
                        # this done for only contain one unique this unique num ber is started by m 
                        pickle.dump(PDB_set_info , open(''.join(["amino_acid_",str(a),"PDB_set_info.p"]), "wb")) 
                        p = p+1
                        checked_group.append(Total_count_group)
                        checked_group.append(p)
                        checked_group.append(staisfied_PDB_ids)
                        #the over all hits are represented in the pikle file name
                        m = m + 1
                        
                        os.chdir('/')
                        os.chdir(self.saving_dir)
                        pickle.dump(checked_group , open(''.join(["amino_acid_",str(a),"_group_",str(m),"_",str(Total_count_group),"_total_hits",".p"]), "wb")) 
                        print("progress unique morethan: ",m)
                        if highest_hit < Total_count_group:
                            highest_hit = Total_count_group
                            #highest_hit_group = []
                            #highest_hit_group.append(m) 
                        if highest_hit == Total_count_group:
                           #highest_hit_group.append(m) 
                           print("m",m) 
                        print("Amino acid checking: ", a, " Coressponding i: ",str(i), " highest_hit achived: ", str(highest_hit))
                        os.chdir('/')
                        os.chdir(self.working_dir)
        print("Amino acid checking: ", a, " Done ")#," highest hit: ", highest_hit, "Groups that has that hits: ", highest_hit_group)            
        #            os.chdir(self.saving_dir)    
        
        
            