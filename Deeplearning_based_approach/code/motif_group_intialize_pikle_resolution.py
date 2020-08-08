# -*- coding: utf-8 -*-
"""
Created on %20-Apr-2018

@author: %A.Nishanth C00294860
"""

import csv
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
import math
import copy#to make deep copy
#%%
"""
896 th line new class started for make property
"""
# to compare the two lists
class motif_group_intialize_pikle:
    """
    This class is basically written for creating the pikle of groups of
    aminoacids with neighbours deatlis with ascending order of length in the group

    """
    def __init__(self,working_dir,saving_dir,pdb_name,resolution):
        
        self.working_dir=working_dir
        self.saving_dir =saving_dir
        self.pdb_name = pdb_name   
        self.resolution = resolution
        
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
    
    #% function for find the equation of triangle and the distnce from the surface
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
    
    #% function for theta angle issue
    def border_angle_issue_solved(self,adding_angle,selected_point_theta_angle,triangle_theta):
        """
        Inputs:
            adding_angle                :  The range we are gonna check
            selected_point_theta_angle  :  Missing point theta angle in radian
            triangle_theta              :  Triangle theta points
            
        Output:
            Atleast one of the theta angle satisified the given criteria(Booleaan)
            
        Inclusion of the fuction:
        Since the radians near to border like -180 degree has to connected with +180 degree(that where the problem caused)
        Inorder to do that only theta angle is canged fro -180-+180 degree
        And Phi angle only changed between -90 to +90 degree(this doesn't have to changed because it decied in top or bottom)
        
        Here everything is in radian
        180 degree = 3.14159
        
        It is imporatnt to make a decider to check the range of the theta angle fell
        """ 
    #    print("work here")
        range_positive = selected_point_theta_angle + adding_angle   
        range_negative = selected_point_theta_angle - adding_angle              
        """
        higher than 3.14 has to changed in to negative range that is 
        range_positive > 3.14
            fell in negative theta range
                -3.14 to  range_positive-3.14
            fell in positive range that is from 
                negative range to 3.14
        
        range_neagtive < -3.14
            fell in negative theta range
                range_positive to -3.14
            fell in positive range that is from negative range to 3.14
                -3.14 to   
        """
        c=0#condtion checker whether that angle satisfied or not
        for chk in triangle_theta:#chk #is the checking angle
            if range_positive > 3.14:
                #fell in negative theta range
                if -3.14 <= chk <= -3.14 + range_positive-3.14:
                    c=1
                #fell in positive range 
                if range_negative <= chk <=3.14:
                     c=1
            elif range_negative < -3.14:
                #fell in negative theta range
                if -3.14 <= chk <= range_positive:
                    c=1
                #fell in positive range 
                # since (range_negative+3.14) is negative number
                if 3.14 +(range_negative+3.14) <= chk <= 3.14:
                     c=1
            elif range_negative <= chk <= range_positive: 
                    c=1
        if c==1:# atleeast once the triangle satisfied any of the criteria
            return True

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
                       adding_angle = 0.2618
                       if self.border_angle_issue_solved(adding_angle,selected_point_theta_angle, theta_t[i]):
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
                       if self.border_angle_issue_solved(0.2618,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(0.5236 ,selected_point_theta_angle, theta_t[i]):
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
                       if self.border_angle_issue_solved(0.2618,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(0.5236 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(0.7853 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                                
            if len(triangle_sat)==0:
                print("PDB-- 120 degree condition: ",self.pdb_name)
                print("corresponding i:::::::::: ",i)
                print("corresponding missing point: ",checking_points_order[s])
                for i in range(0,len(phi_t)):
                    phi_condition=0
                    for o in range(0,3):
                        if selected_point_phi_angle - 1.0472 <= phi_t[i][o] <= selected_point_phi_angle + 1.0472  :
                            phi_condition=1
                    if phi_condition == 1:
                       if self.border_angle_issue_solved(0.2618,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(0.5236 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(0.7853 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(1.0472 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                                
            if len(triangle_sat)==0:
                print("PDB-- 150 degree condition: ",self.pdb_name)
                print("corresponding i:::::::::: ",i)
                print("corresponding missing point: ",checking_points_order[s])
                for i in range(0,len(phi_t)):
                    phi_condition=0
                    for o in range(0,3):
                        if selected_point_phi_angle - 1.309 <= phi_t[i][o] <= selected_point_phi_angle + 1.309  :
                            phi_condition=1
                    if phi_condition == 1:
                       if self.border_angle_issue_solved(0.2618,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(0.5236 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(0.7853 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(1.0472 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(1.309 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                                
                       
            if len(triangle_sat)==0:
                print("PDB-- 180 degree condition: ",self.pdb_name)
                print("corresponding i:::::::::: ",i)
                print("corresponding missing point: ",checking_points_order[s])
                for i in range(0,len(phi_t)):
                    phi_condition=0
                    for o in range(0,3):
                        if selected_point_phi_angle - 1.5706  <= phi_t[i][o] <= selected_point_phi_angle + 1.5706:
                            phi_condition=1
#                            print("Phi condition satisfied")
                    if phi_condition == 1:
                       if self.border_angle_issue_solved(0.2618,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(0.5236 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(0.7853 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(1.0472 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                       elif self.border_angle_issue_solved(1.5706 ,selected_point_theta_angle, theta_t[i]):
                            triangle_sat.append(i)
                                
            if len(triangle_sat)==0:
                print("-------Only phi condition chk----: ",self.pdb_name)
                print("corresponding i:::::::::: ",i)
                print("corresponding missing point: ",checking_points_order[s])
                raise ValueError("Problem encounterd in surface triangle formation")

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
            
            since these triangle_break and theta and phi angles are used again and again after updated also
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
        self.coordinates = coordinates/(self.resolution)     
        self.coordinates_polar = coordinates_polar      
        self.tri = tri
        self.vertices = verrtices
#        return  coordinates, coordinates_polar, tri, verrtices  
        
    def distance_edge(self,uncovered_point,un_covered_edges):
        """
        This function calculates the perpendicular diatance from the uncovered point to the uncovered edge
        And return the nearest edge
        
        uncovered_point : The point wanted to find the nearest edge
        un_covered_edges: Uncovered edges
        """
        coordinates  = self.coordinates  
        missed_point = coordinates[uncovered_point]   
        # go through the uncovered edges
        distances =[]
        for i in range(0,len(un_covered_edges)):
            p1 = coordinates[un_covered_edges[i][0]]
            p2 = coordinates[un_covered_edges[i][1]]
            d = np.linalg.norm(np.cross(p2-p1, p1- missed_point))/np.linalg.norm(p2-p1)
            if d<0:
                d=(-1)*d
            distances.append(d)
        distances.sort()
        for i in range(0,len(un_covered_edges)):
            p1 = coordinates[un_covered_edges[i][0]]
            p2 = coordinates[un_covered_edges[i][1]]
            d = np.linalg.norm(np.cross(p2-p1, p1- missed_point))/np.linalg.norm(p2-p1)
            if d<0:
                d=(-1)*d
            if distances[0]==d:# check the minimum distance edge
               return copy.deepcopy(un_covered_edges[i])
    
    def edge_length(self, p1, p2):
        """
        This one calculate the length between the two points
        """
        return math.sqrt( ( (p2[0] - p1[0]) **2) + ( (p2[1] - p1[1]) **2) +( (p2[2] - p1[2]) **2)  ) 
    
    
    
    def dip_area_volume(self, tri_t, dip_chosen, coordinates):
        """
        Function calulate the area of the triangle
        Atoms of the surface number
        coordinates of all surface atoms
        
        First calculate the edge length of the triangles seperately and save them as a, b, c
        
        Then from that triangle find the volume of the tetrahydron formed 
        """
        a = self.edge_length(coordinates[tri_t[0]], coordinates[tri_t[1]])
        b = self.edge_length(coordinates[tri_t[0]], coordinates[tri_t[2]])    
        c = self.edge_length(coordinates[tri_t[1]], coordinates[tri_t[2]])   
    
        # calculate the semi-perimeter
        s = (a + b + c) / 2
        area = math.sqrt((s*(s-a)*(s-b)*(s-c)))
       
        d = self.distance_surface(coordinates,tri_t,dip_chosen)#for finding the distance from the surface
        volume = (1/3)*area*d
        return area,volume
    
    def triangle_uncovered_edge(self,un_covered_edges,tri_t):
        """
        Change the given triangle into uncovered edges
        And each of these edges are in the ascending order
        """
        a = [tri_t[0],tri_t[1]]
        b = [tri_t[0],tri_t[2]]
        c = [tri_t[1],tri_t[2]]
        a.sort()
        b.sort()
        c.sort()
        un_covered_edges.append(a)
        un_covered_edges.append(b)
        un_covered_edges.append(c)
        return un_covered_edges
    
    def center_atom_group(self,center_atom_index,group_from_coordinates):
        """
        It calculates the atoms formed as group by the center atom
        center_atom_index: index of the center atom
        group_from_coordinates: conatin other traingle vertices with out the center atom
        """
        atoms_in_group = []
        for i in  range(0,len(group_from_coordinates[center_atom_index])):
            atoms_in_group.append(group_from_coordinates[center_atom_index][i][0])
            atoms_in_group.append(group_from_coordinates[center_atom_index][i][1])
        atoms_in_group = list(set(atoms_in_group))
        return atoms_in_group 
    
    def choose_dip_bump_struct(self):
        """
        It find out the substructure which is dip or dump

        Inorder to do this group_from_coordinates is used
        group_from_coordinates each surface atom is once considerd as center atom and their corresponding substructures are formed from that
        these groups are reprecented by the index of the surface atoms
        
        Calculate the normalized coordinate by divide the coordinates by resolution
        -calculate the center of gravity from the newly calculated the coordinates
        
        choose the center atom once
        -first find the center atom is: dip or bump
            where   dip means the radius of center is less than all the other radius from the gravity point
                    bump means the radius of center is higher than all the other radius from the gravity point
        
        these dip and dump data are saved differently to train the model differently
        

        ------------------------------------------------------------------------------------------------------------
        
        Bump case:
            Surface:
                Just take the surface of triangles with the center
                and add them
                
            Volume:    
                Order the radius of the sorrunding atoms in assending order
                choose the triangles in that order and form  tetrahedron with the center atom
                add them up
                
            Depth:
        Perpendicular distance from the first trangle found for surface
            """
        
        coordinates_polar = self.coordinates_polar
        coordinates = self.coordinates
        verrtices = self.vertices
        """
        go through the surface atoms and findoutthe hit map of triangles in triangles
        """
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
        
        
        #% choose the center atom once
        dip_center_atoms =[]
        bump_center_atoms =[] 
        #-first find the center atom is: dip or bump
        # inorder to do thi suse the polar coordinate data radius
        for center_atom_index in range(0,len(coordinates_polar)):
            center_radius = coordinates_polar[center_atom_index][0]
            atoms_in_group_c = self.center_atom_group(center_atom_index,group_from_coordinates)
           
            # if these conditions are in one means these condition staisfied
            # if these conditions changed to 0 means it broken
            dip_cond = 1
            bump_cond = 1
            for i in range(0,len(atoms_in_group_c)):
                if center_radius < coordinates_polar[i][0]:
                    bump_cond = 0#bump condition broken
                elif center_radius > coordinates_polar[i][0]:
                    dip_cond  = 0#dip condition broken
                    
            if dip_cond == 1:
                dip_center_atoms.append(center_atom_index)
            elif bump_cond == 1:
                bump_center_atoms.append(center_atom_index)
        
        self.dip_center_atoms = dip_center_atoms
        self.bump_center_atoms = bump_center_atoms
        self.group_from_coordinates =  group_from_coordinates
        self.amino_acid_details = self.aminoacid_function_retrieve(''.join([self.pdb_name,'_aminoacids.csv']))
        
    def dip_group_amino_acid_pikle(self):
        """
        It includes the details of the area and volume 
        
        Dip case:
            surface:
                first take the sorrounding atoms(other than center atom) order them in desending order with radius
                then make triangles from the order and consider previously made traiangle(edges)
                Then add all the triangle surfaces.
                
            Volume:    
                find the volume for each tetrahedron(traingle formed for surface with center atom)
                add those volumes
                
            Depth:
                Perpendicular distance from the first trangle found for surface
                
        """
        
        dip_center_atoms = self.dip_center_atoms 
        group_from_coordinates = self.group_from_coordinates
        coordinates_polar = self.coordinates_polar
        coordinates = self.coordinates
        amino_acid_details = self.amino_acid_details
        hashmap_edge = self.hashmap_edge
        
        group = []# to store the frop details PDB file
        group_details=[]
#        print("Here")
        for dip_chosen in dip_center_atoms:
            #  first take the sorrounding atoms(other than center atom) order them in desending order with radius
            atoms_in_group_c = self.center_atom_group(dip_chosen,group_from_coordinates)
            atoms_in_group_radius = []
            for i in range(0,len(atoms_in_group_c)):
                atoms_in_group_radius.append(coordinates_polar[atoms_in_group_c[i]][0])
            
            atoms_in_group_radius.sort(reverse=True)  
            atoms_in_group_in_radius_order =[]
            for i in atoms_in_group_radius:
                for j in range(0,len(atoms_in_group_c)):
                    if i == coordinates_polar[atoms_in_group_c[j]][0]:
                        atoms_in_group_in_radius_order.append(atoms_in_group_c[j])
                        
            # then make triangles from the order and consider previously made traiangle(edges)
            triangles_surface = []
            tri_t =atoms_in_group_in_radius_order[0:3]
            triangles_surface.append(tri_t)
            
            """
            Inorder to calculate the depth 
            
            Center atom dip chosen to the highest elevation triangle surface is calculated
            """
            depth = self.distance_surface(coordinates,tri_t,dip_chosen)
#            print("depth")
#            print("depth: ",depth)
            uncovered_points=[]
            if len(atoms_in_group_in_radius_order) > 3:
                uncovered_points = atoms_in_group_in_radius_order[3:len(atoms_in_group_in_radius_order)]
            
            #% then take that triangle and choose the next triangle
            """
            Inorder to make trangulation of surface atoms
            -Findout the uncovered atoms in the group
            -if uncovered exist 
            Take the uncovered triangle edges one by one and check the uncovered point near to which
             edge(perpendicular distance from the edge) and assign that as 
            triangle and remove the edge from uncovered edge and the point from uncovered point
            
            """
            un_covered_edges = []
            un_covered_edges =  self.triangle_uncovered_edge(un_covered_edges,tri_t)
            
            #  Then add all the triangle surfaces.
            while len(uncovered_points)>0:
                uncovered_point = uncovered_points[0]
                near_edge = self.distance_edge(uncovered_point,un_covered_edges)
                #then form a triangle from that and attach it
                tri_t = copy.deepcopy(near_edge)
                tri_t.append(uncovered_point)
                triangles_surface.append(tri_t)
                un_covered_edges =  self.triangle_uncovered_edge(un_covered_edges,tri_t)
                #then remove  the nearest edge and the uncovered edge from the list
                # because it occurs earlir and now
                un_covered_edges.remove(near_edge)
                un_covered_edges.remove(near_edge)
                #% remove the uncovered point from the uncovered_points
                uncovered_points.remove(uncovered_point)
            #print("Done: ",triangles_surface)
            
            #% calculate the area and volume of the substructure
            #% then find the volume and the surface from thge triangles
            
            total_area = 0
            total_volume= 0
            for tri_t in triangles_surface:
                area,volume =self.dip_area_volume(tri_t, dip_chosen, coordinates)
                total_area = total_area + area
                total_volume = total_volume + volume
            
            #% then change the values to aminoacids
            temp=[]
            # first append the center atom amino acid details
            for chk_group in group_from_coordinates[dip_chosen]:
                #take the group detail and find out the index of aminoacid
                index_1 = self.amino_acid_val(amino_acid_details[chk_group[0]])
                index_2 = self.amino_acid_val(amino_acid_details[chk_group[1]])
                temp.append(self.hashmap_edge[index_1][index_2])
            temp.sort()# because this only have the neighbour thus easy to findout the same group
            temp_2 = []
            temp_2.append(hashmap_edge[self.amino_acid_val(amino_acid_details[dip_chosen])][self.amino_acid_val(amino_acid_details[dip_chosen])])
            temp_2.append([depth, total_volume, total_area])
            temp_2.append(temp)
            group.append(temp_2)
            group_details.append([self.amino_acid_val(amino_acid_details[dip_chosen]),len(temp)])
            
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
            #to avoid save redundant senter amino acid group    
            if len(sub_group_temp)>0:
                os.chdir('/')
                os.chdir(self.saving_dir)
    #            pickle.dump(group, open(''.join([self.pdb_name,"_dip_group_",str(j),".p"]), "wb"))
    #            pickle.dump(group_details, open(''.join([self.pdb_name,"_dip_group_details",str(j),".p"]), "wb"))
                pickle.dump(sub_group_temp, open(''.join([self.pdb_name,"_dip_group_",str(j),".p"]), "wb")) 
            
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


class motif_group_pikle_to_property_pikle:
    """
    This class is basically written for reading the pikle of groups and make the property values for that PDB
    """
    def __init__(self, working_dir, saving_dir, group_files, pdb_name_unique, pdb_name):
        
        self.working_dir=working_dir
        self.saving_dir =saving_dir
        self.pdb_name = pdb_name   
        '''
        Since this is loading the previous done pikle files 
        it need to has the list pdb file should be done earlier
        '''
        self.group_files = group_files
         
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

    #% assigning property values
    def amino_property(self,index):
        """
        amino acid property values from their indexes
        """
        if (index ==0):
            amini_acid_properties=np.array([0,0,1,1,0,0,0,0,0,0,0,0,1,0,0.404040404,0.087621697])
        elif (index ==1):
            amini_acid_properties=np.array([0,1,0,0,1,0,1,0,0,0,0,0,1,0,0.207070707,0.079276773])
        elif(index ==2):
            amini_acid_properties=np.array([1,0,0,0,0,1,0,0,0,0,1,0,0,1,0.146464646,0.065368567])
        elif(index ==3):
            amini_acid_properties=np.array([1,0,0,0,0,1,0,0,0,0,1,0,0,1,0.747474747,1])
        elif(index ==4):
            amini_acid_properties=np.array([1,0,0,0,0,1,0,0,0,1,0,1,0,0,0.404040404,0.02364395])
        elif(index ==5):
            amini_acid_properties=np.array([1,0,0,0,0,1,0,0,0,1,0,1,0,0,0.439393939,0.066759388])
        elif(index ==6):
            amini_acid_properties=np.array([0,1,0,0,0,1,1,0,0,0,0,0,1,0,0.166666667,0.063977747])
        elif(index ==7):
            amini_acid_properties=np.array([0,1,0,0,0,1,1,0,0,0,0,0,1,0,0,0.043115438])
        elif(index ==8):
            amini_acid_properties=np.array([0,1,0,0,1,0,0,0,0,0,1,0,0,1,0.085858586,0.009735744])
        elif(index ==9):
            amini_acid_properties=np.array([0,1,0,0,0,1,1,0,0,0,0,0,1,0,0.176767677,0.069541029])
        elif(index ==10):
            amini_acid_properties=np.array([0,1,0,0,0,1,1,0,0,0,0,0,1,0,0.161616162,0.061196106])
        elif(index ==11):
            amini_acid_properties=np.array([0,1,0,1,0,0,0,1,0,0,0,0,1,0,0.156565657,0.068150209])
        elif(index ==12):
            amini_acid_properties=np.array([0,1,0,0,1,0,1,0,0,0,0,0,1,0,1,0])
        elif(index ==13):
            amini_acid_properties=np.array([0,1,0,1,0,0,0,1,0,0,0,0,1,0,0.297979798,0.093184979])
        elif(index ==14):
            amini_acid_properties=np.array([0,0,1,1,0,0,0,0,1,0,0,0,1,0,0.54040404,0.089012517])
        elif(index ==15):
            amini_acid_properties=np.array([0,0,1,1,0,0,0,0,1,0,0,0,1,0,0.484848485,0.084840056])
        elif(index ==16):
            amini_acid_properties=np.array([0,0,1,1,0,0,0,0,1,0,0,0,1,0,0.404040404,0.090403338])
        elif(index ==17):
            amini_acid_properties=np.array([0,0,1,1,0,0,0,1,0,0,0,0,1,0,0.222222222,0.121001391])
        elif(index ==18):
            amini_acid_properties=np.array([0,0,1,1,0,0,0,0,1,0,0,0,1,0,0.464646465,0.080667594])
        elif(index ==19):
            amini_acid_properties=np.array([0,0,1,1,0,0,0,0,0,0,0,0,1,0,0.909090909,0.038942976])
            
        return amini_acid_properties      
    
    #% then add the properties of the sub_group
    def intialize_sub_group_property_matrix(self):
        """
        Initialization of 20 amino acid properties as metrix inorder to access later
        This sub_group_properties_intialized can be used for all subgroups because it remain same
        """
        # create an empty matrix to place the properties
        sub_group_properties_intialized = np.zeros((20,16))           
        # asssign the properties to all
        for i in range(0,len(sub_group_properties_intialized)):
            sub_group_properties_intialized[i][:] = self.amino_property(i)
        self.sub_group_properties_intialized = sub_group_properties_intialized



    #% hash map to index
    def hash_to_indexes(self,chk_point):
        """
        This check the chk_point value belongs to which indices
        """
        for i in range(0,20):
          for j in range(0,20):
            if j-i-1<0:
                if chk_point == self.hashmap_edge[i,j]:
                    return [i,j]

     
    """
    This extract the information from the given pikle file
    These info can be used in different ways reprecentation
    From here extract the general information by the functions
            
    since these pikles are varying in size so functions use use list to store their values
    
    These functions
        pikle_inform_decode
        find_vertices_in_group
                                 can be used for other purposes
    """
    def pikle_inform_decode(self,pikle_name):
        """
        Retrieve the information in the format from the pikle file 
        
        group_indexes_t: this save the subgroup format into
            [depth, total_volume, total_area]
            center_atom_property_index 
            group_of_edge_index like index_1 and index_2
            
        group_indexes_temp: contain all sub_groups(group_indexes_t) of one pikle group
        """
        os.chdir('/')
        os.chdir(self.working_dir)
        group_indexes_temp =[]
        group_chk= pickle.load(open(pikle_name, "rb"))
        center_atom = group_chk[0][0]
        for i in range(0,len(group_chk)):
            group_indexes_t=[]
            group_indexes_t.append(group_chk[i][1])
            group_indexes_t.append(center_atom)
            for chk_group in group_chk[i][2]:
                group_indexes_t.append(self.hash_to_indexes(chk_group))  
            group_indexes_temp.append(group_indexes_t)
        return group_indexes_temp
    
    def find_vertices_in_group(self,sub_group_temp):
        """
        Find out the vertices in the group
        This function reduce the information content 
        Because when it reduce to edges to vertices 
        the information about the vertice connection lossed
        
        sub_group_temp: Contain the edges
        return:
            vertices: are the idices of amino acids for the property
            count: number of time occurs in the subgroup 
                    respect to vertices
        """
        #% First find out the emlements on the group
        group_elements=[]
        for edge in sub_group_temp:
            group_elements.append(edge[0])
            group_elements.append(edge[1])
            
        unique_sub_group =  list(set(group_elements))
        unique_sub_group_count = []
        for i in range(0,len(unique_sub_group)):
            c=0
            for ele in  group_elements:
                if unique_sub_group[i]==ele:
                    c=c+1
            #c devided by two because one vertice has two edge connection
            unique_sub_group_count.append(int(c/2))     
        return unique_sub_group,unique_sub_group_count    

    def group_temp_properties(self,sub_group_temp):
        """
        It doesn't consider center atom property
        Function create a 20 x 16 property matrix
        Where the 
            each row present the aminoacid
            coloum represent the property
        Input:
            sub_group_temp: contain the edges in the sub_group without the center atom
            
        o/p: 20 x 16 
            each group multiplied by the factor of howmany time it occured in the group
        """
        sub_group_properties_intialized = self.sub_group_properties_intialized
        
        # first intialize the subgroup to contain aminoacid properties
        sub_group_temp_properties = copy.deepcopy(sub_group_properties_intialized) 
        #% first make the index multiplying values
        unique_sub_group,unique_sub_group_count = self.find_vertices_in_group(sub_group_temp)
        # then create a multiplying factor
        multiply_factor = np.zeros((1,len(sub_group_properties_intialized)))
        for i in range(0,len(unique_sub_group)):
            multiply_factor[0][unique_sub_group[i]] = unique_sub_group_count[i]
        # then multiply the factor to the properties
        for i in range(0,len(sub_group_properties_intialized)):
            sub_group_temp_properties[i][:] = np.dot(sub_group_temp_properties[i][:],multiply_factor[0][i])
        #    print(multiply_factor[0][i])
        return sub_group_temp_properties
 
    """
    This is the main function which make property value
    This way of representation can be changed
    """
    def property_flatten_array(self,d_v_a, center_amino_acid_property, sub_group_temp_properties):
        """
        Input: 
            d_v_a                       : Depth, Volume, and Area information in order as numpy array
            center_amino_acid_property  : Amino acid in the center dip's property values
            sub_group_temp_properties   : Ami no acid properties in dip with out center amino acid properties
        
        o/p:
            Combine the properties with depth, Volume and, Area
            
        Make the property matrix combined with 
            area, volume, Depth
            Here the center atom, and the sorrunding atoms are treated in different way
            
        The properties are attached in the order as follows
            The depth, volume, area are attached
                3
                
            Center atoms property is multiplied by 
                (1,depth, 1/depth, Volume, 1/Volume, 1/Area, Area)
                1 x 16 x 7 = 112
                
            Sorrunding amino acid properties  are make like
                multipled by
                (1, Volume, 1/Volume, 1/Area, Area)
                20 x 16 x 5 = 1600
                        
            So total features 
                1600 + 112 + 3 = 1715
                Less than (42 x 42 image pixel feature with one channel)
                          (24 x 24 image pixel feature with three channels)
                          
                          MNIST: 784(28 x 28 with one channel)
                          
            Case haven't considered for time constraint:
                Area or 
                Volume or 
                depth 
                    less than one is treated the same way as values higher than one
        """
        
        sub_group_properties = np.zeros((1,1715))
        sub_group_properties[0][0:3] = copy.deepcopy(d_v_a)
        sub_group_properties[0][3:19] = copy.deepcopy(center_amino_acid_property)
        
        """
        Center atoms property is multiplied by 
                (1,depth, 1/depth, Volume, 1/Volume, 1/Area, Area)
        """
        sub_group_properties[0][19:35] = np.dot(sub_group_properties[0][3:19],sub_group_properties[0][0])
        sub_group_properties[0][35:51] = np.dot(sub_group_properties[0][3:19],1/sub_group_properties[0][0])
        sub_group_properties[0][51:67] = np.dot(sub_group_properties[0][3:19],sub_group_properties[0][1])
        sub_group_properties[0][67:83] = np.dot(sub_group_properties[0][3:19],1/sub_group_properties[0][1])
        sub_group_properties[0][83:99] = np.dot(sub_group_properties[0][3:19],sub_group_properties[0][2])
        sub_group_properties[0][99:115] = np.dot(sub_group_properties[0][3:19],1/sub_group_properties[0][2])
        """
        Sorrunding amino acid properties  are make like
                multipled by
                (1, Volume, 1/Volume, 1/Area, Area)
        """
        sub_group_temp_property_flatten = copy.deepcopy(sub_group_temp_properties.flatten())
        sub_group_properties[0][115:435]= sub_group_temp_property_flatten
        sub_group_properties[0][435:755]= np.dot(sub_group_properties[0][115:435],sub_group_properties[0][1])
        sub_group_properties[0][755:1075]= np.dot(sub_group_properties[0][115:435],1/sub_group_properties[0][1])
        sub_group_properties[0][1075:1395]= np.dot(sub_group_properties[0][115:435],sub_group_properties[0][2])
        sub_group_properties[0][1395:1715]= np.dot(sub_group_properties[0][115:435],1/sub_group_properties[0][2])
        
        return sub_group_properties



    def main_fn_property_from_PDB(self):
        """
        Main function call other helper functions tpo make property files as pikle
        """
        #first findout the which pikle files belong to the given PDB file
        chk_pdb_group =[]
        for l in self.group_files:
            if self.pdb_name  in l:
                chk_pdb_group.append(self.pikle_inform_decode(l))  
        """
        Initial format of the property file
            Each group has 
                
                to represent the aminoacids in the group 20 x 16 
            Then take each PDB ids groups sepeartely and calculate their
            Property values
            
            The property of the values are represented by 20 x 16 which tells serrounded amino acid properies added

        """       
                        
        #% using the  group_indexes_temp create the pikle file of the property
        property_all_pdb =[]
        # then choose one of the center atoms pilkle group
        for i in range(0,len(chk_pdb_group)):
            #since all in this group center amino acid property is same
            center_amino_acid_property = self.amino_property(self.hash_to_indexes(chk_pdb_group[i][0][1])[0])
            for j in range(0,len(chk_pdb_group[i])):
                sub_group_temp = []
                ''' 
                depth, volume, area information
                '''
                d_v_a = np.array([chk_pdb_group[i][j][0]])
                for k in range(2,len(chk_pdb_group[i][j])):
                    sub_group_temp.append(chk_pdb_group[i][j][k])
                sub_group_temp_properties = self.group_temp_properties(sub_group_temp)
                property_all_pdb.append(self.property_flatten_array(d_v_a, center_amino_acid_property, sub_group_temp_properties))
        
        """
        Save as pikle_file
        """
        os.chdir('/')
        os.chdir(self.saving_dir)
        pickle.dump(property_all_pdb , open(''.join(["property_all_pdb_",self.pdb_name,".p"]), "wb")) 
        print("progress pdb id: ",self.pdb_name)
        print("")
        os.chdir('/')
        os.chdir(self.working_dir)
        
        
class motif_group_find_pikle:
    """
    But not editted with resolution dip results checked
    
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
                            highest_hit_group = []
                            highest_hit_group.append(m)
                        if highest_hit == Total_count_group:
                           highest_hit_group.append(m) 
                        print("Amino acid checking: ", a, " Coressponding i: ",str(i), " highest_hit achived: ", str(highest_hit))
                        os.chdir('/')
                        os.chdir(self.working_dir)
        if highest_hit == 0:
            print("Amino acid checking: ",a," Nothing_satisfied_this_condition")
            print("---------------------------Done----------------------------")
        else:
            print("Amino acid checking: ", a, " Done .... ----- ddddddddooooooooooonnnnnnnneeeeeeee"," hoighest hit: ", highest_hit, "Groups that has that hits: ", highest_hit_group)            
        #            os.chdir(self.saving_dir)    
        
        
            
