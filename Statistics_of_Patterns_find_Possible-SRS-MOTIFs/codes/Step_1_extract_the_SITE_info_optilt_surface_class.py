# -*- coding: utf-8 -*-
"""
Created on %%(26-June-2020) at 11.09 A.m

@author: %A.Nishanth C00294860
corrected on 11-Oct-2020 at 5.28p.m
"""

import os
import pickle
from copy import deepcopy
import numpy as np

class surface_msms_depth_MOTIF__extract_class_quat:
    """
    USE MSMS tool(1996) to calculte the soluable access area of the surface
        depth of the C_alpha carbons, and 
              of the residue 
                          from the surface
    this vertion fix the principle direction
    """
    def __init__(self, loading_pikle_dir, saving_dir, pdb_name,surf_threh_cond=True):
        self.saving_dir = saving_dir
        self.pdb_name  = pdb_name #name of the class
       
        '''Assign property way'''
        self.old_17_prop=True
        #to select the optimal direction it should be true
        self.optimal_tilt=True
        '''loading_the_thresholds
        Earlier thresh_hold_ca=7.2,thresh_hold_res=6.7
        '''
#        os.chdir("/")
#        os.chdir(loading_dir_threshold)
        #        thresh_hold_ca= pickle.load(open("max_depths_ca_MOTIF.p", "rb"))  
        #        thresh_hold_res=pickle.load(open("max_depths_res_MOTIF.p", "rb"))  
        thresh_hold_ca=7.2
        thresh_hold_res=6.7

        os.chdir('/')
        os.chdir(loading_pikle_dir)
        coordinates = pickle.load( open( ''.join(["coordinates_",pdb_name,".p"]), "rb" ))
        aminoacids = pickle.load( open( ''.join(["amino_acid_",pdb_name,".p"]), "rb" ))
        fin_res_depth_all = pickle.load( open( ''.join(["fin_res_depth_all_",pdb_name,".p"]), "rb" ))
        fin_ca_depth_all = pickle.load( open( ''.join(["fin_ca_depth_all_",pdb_name,".p"]), "rb" ))
        MOTIF_indexs_all = pickle.load( open( ''.join(["MOTIF_indexs_all_",pdb_name,".p"]), "rb" ))         

        c_alpha_indexes_MOTIF= sum(MOTIF_indexs_all, [])
        res_factor = 2.25 # see the documentation twhy 2.25 is chosen

        sur_res = []
        sur_res_cor_intial = []
        MOTIF_prop =[]
        MOTIF_prop_t=[]
        #% to find out the surface atoms residues 
        for i in range(0,len(fin_res_depth_all)):
            if surf_threh_cond:
                if fin_ca_depth_all[i] <= thresh_hold_ca:
                    if fin_res_depth_all[i] <= thresh_hold_res:
                        sur_res.append(aminoacids[i])
                        # multiply each coordinate by 2 (just for increasing the resolution) and then round them to decimal numbers.
                        #sur_res_cor_intial_round.append([round(res_factor*coordinates[i][0]),round(res_factor*coordinates[i][1]),round(res_factor*coordinates[i][2])])
                        sur_res_cor_intial.append([res_factor*coordinates[i][0],res_factor*coordinates[i][1],res_factor*coordinates[i][2]])
    
                        if i in c_alpha_indexes_MOTIF:
                            for j in range(0,len(MOTIF_indexs_all)):
                                if i in MOTIF_indexs_all[j]:
                                    MOTIF_prop.append(1+j)
                                    MOTIF_prop_t.append(1)
                        else:
                            MOTIF_prop.append(0)
                            MOTIF_prop_t.append(0)

            else:
                sur_res.append(aminoacids[i])
                # multiply each coordinate by 2 (just for increasing the resolution) and then round them to decimal numbers.
                #sur_res_cor_intial_round.append([round(res_factor*coordinates[i][0]),round(res_factor*coordinates[i][1]),round(res_factor*coordinates[i][2])])
                sur_res_cor_intial.append([res_factor*coordinates[i][0],res_factor*coordinates[i][1],res_factor*coordinates[i][2]])

                if i in c_alpha_indexes_MOTIF:
                    for j in range(0,len(MOTIF_indexs_all)):
                        if i in MOTIF_indexs_all[j]:
                            MOTIF_prop.append(1+j)
                            MOTIF_prop_t.append(1)
                else:
                    MOTIF_prop.append(0)
                    MOTIF_prop_t.append(0)

#        if len(c_alpha_indexes_MOTIF)!=np.sum(MOTIF_prop):
#            print(len(c_alpha_indexes_MOTIF)-np.sum(MOTIF_prop)," MOTIF atoms miss in ",pdb_name)
        if len(sur_res_cor_intial)==0:       
            print(pdb_name, ' not has single atom to satify')
        else:
#            '''first title the coordinates in principle direction; and finally round the o/p coordinates'''
#            sur_res_cor_tilted = self.change_C_alpha_principle_direc(sur_res_cor_intial)
            '''rotate in optimal direction'''
            sur_res_cor_tilted = self.rotate_C_alpha_optimal_direc(sur_res_cor_intial)

            '''
            just take the optimally tilted coordinates 
            And extract only the MOTIF coordinates information
            '''
            MOTIF_tilted_coordinates=np.empty((np.sum(MOTIF_prop_t),3))
            MOTIF_res=[]
            MOTIF_group_info=[]
            j=0
            for i in range(0,len(sur_res_cor_tilted)):
                if MOTIF_prop[i]!=0:
                    MOTIF_tilted_coordinates[j,:] = sur_res_cor_tilted[i]
                    MOTIF_res.append(sur_res[i])
                    MOTIF_group_info.append(MOTIF_prop[i])
                    j=j+1
            _,MOTIF_quarter = self.quaterize_coordinate(MOTIF_tilted_coordinates)
           
            os.chdir('/')
            os.chdir(saving_dir)
            pickle.dump(MOTIF_tilted_coordinates, open( ''.join(["MOTIF_tilted_coordinates_",self.pdb_name,".p"]), "wb" ) ) 
            pickle.dump(MOTIF_res, open( ''.join(["MOTIF_res_",self.pdb_name,".p"]), "wb" ) ) 
            pickle.dump(MOTIF_quarter, open( ''.join(["MOTIF_quarter_",self.pdb_name,".p"]), "wb" ) ) 
            pickle.dump(MOTIF_group_info, open( ''.join(["MOTIF_group_info_",self.pdb_name,".p"]), "wb" ) ) 
 
    def quaterize_coordinate(self,coordinates):
        """ 
        Take the coordinates and find which quarter it fell
        
        find the minimum position to move this shift the coordinates to make NN-train
        If the coordinate fell in quarter 1, assign 1 in quarter 1 for that coordinate
        
        zero_cl=give some gap to the quarters overlap to use the continuity between them
        
        prefered -5 to use if the kernal size 7 is used 
        """
        zero_cl=-5
        quarter=np.zeros((8,len(coordinates)))

        s_coordinates = deepcopy(coordinates)
        
        for i in range(0,len(coordinates)):
            if coordinates[i][2] >= zero_cl:
                s_coordinates[i][2] = coordinates[i][2]-zero_cl

                if coordinates[i][0] >= zero_cl and coordinates[i][1] >=zero_cl:
                    quarter[0][i]=1
                    s_coordinates[i][0] = coordinates[i][0]-zero_cl
                    s_coordinates[i][1] = coordinates[i][1]-zero_cl                 
                elif coordinates[i][0] < -zero_cl and coordinates[i][1] >= zero_cl:
                    quarter[1][i]=1
                    s_coordinates[i][0] = (-1)*(coordinates[i][0]+zero_cl)
                    s_coordinates[i][1] = coordinates[i][1]-zero_cl                                        
                elif coordinates[i][0] < -zero_cl and coordinates[i][1] < -zero_cl:
                    quarter[2][i]=1
                    s_coordinates[i][0] = (-1)*(coordinates[i][0]+zero_cl)
                    s_coordinates[i][1] = (-1)*(coordinates[i][1]+zero_cl)
                elif coordinates[i][0] >= zero_cl and coordinates[i][1] < -zero_cl:
                    quarter[3][i]=1
                    s_coordinates[i][0] = coordinates[i][0]-zero_cl
                    s_coordinates[i][1] = (-1)*(coordinates[i][1]+zero_cl)
                
            elif coordinates[i][2] <  -zero_cl:
                s_coordinates[i][2] = (-1)*(coordinates[i][2]+zero_cl)

                if coordinates[i][0]>=zero_cl and coordinates[i][1] >=zero_cl:
                    quarter[4][i]=1
                    s_coordinates[i][0] = coordinates[i][0]-zero_cl
                    s_coordinates[i][1] = coordinates[i][1]-zero_cl
                elif coordinates[i][0] < -zero_cl and coordinates[i][1] >= zero_cl:
                    quarter[5][i]=1
                    s_coordinates[i][0] = (-1)*(coordinates[i][0]+zero_cl)
                    s_coordinates[i][1] = coordinates[i][1]-zero_cl
                elif coordinates[i][0] < -zero_cl and coordinates[i][1] < -zero_cl:
                    quarter[6][i]=1
                    s_coordinates[i][0] = (-1)*(coordinates[i][0]+zero_cl)
                    s_coordinates[i][1] = (-1)*(coordinates[i][1]+zero_cl)   
                elif coordinates[i][0] >= zero_cl and coordinates[i][1] < -zero_cl:
                    quarter[7][i]=1
                    s_coordinates[i][0] = coordinates[i][0]-zero_cl
                    s_coordinates[i][1] = (-1)*(coordinates[i][1]+zero_cl)
        
        return s_coordinates,quarter
    
  

#    def results(self):
#        '''For creation and checking purpose'''
#        property_surface = self.property_surface
#        MOTIF_prop = self.MOTIF_prop
#        sur_res_cor = self.sur_res_cor
#        return property_surface,MOTIF_prop,sur_res_cor
    
    def shift_coordinate(self,coordinates):
        """ 
        find the minimum position to move this shift the coordinates to make NN-train
        """
        s_coordinates = deepcopy(coordinates)
        m_x = coordinates[0][0]
        m_y = coordinates[0][1] 
        m_z = coordinates[0][2]     
        for i in range(0,len(coordinates)):
            if coordinates[i][0] < m_x:
                m_x = coordinates[i][0]
            if coordinates[i][1] < m_y:
                m_y = coordinates[i][1]
            if coordinates[i][2] < m_z:
                m_z = coordinates[i][2]
        # this will linearly change the position
        for i in range(0,len(coordinates)):
             s_coordinates[i][0]  = coordinates[i][0] - m_x 
             s_coordinates[i][1]  = coordinates[i][1] - m_y
             s_coordinates[i][2]  = coordinates[i][2] - m_z
        return s_coordinates

    """ using eigen vector to rotate in optimal direction"""
    def rotate_C_alpha_optimal_direc(self,coordinates):
        """
        then change the C-alpha coordinates to the optimal direction
        Inorder to do that calculate the eigen vector and 
        use the eigen vector to rotate the coordinates
        
        returns 
        coordinates             : cetralised coordinates
        finalised_cartesian     : Get the optimal direction rotated coordinates
        """
        cen_coordinates = deepcopy(coordinates)
        
        g_x = 0
        g_y = 0
        g_z = 0
        
        x =[]
        y=[]
        z=[]
        for i in range(0,len(coordinates)):
            #find the center of gravity
             g_x = g_x + coordinates[i][0]       
             g_y = g_y + coordinates[i][1]
             g_z = g_z + coordinates[i][2]    
             
             x.append(coordinates[i][0])
             y.append(coordinates[i][1])
             z.append(coordinates[i][2])
        
        #% then centralize the coordinates
        for i in range(0,len(coordinates)):
            #find the center of gravity
             cen_coordinates[i][0]  = coordinates[i][0] - g_x/len(coordinates)    
             cen_coordinates[i][1]  = coordinates[i][1] - g_y/len(coordinates)    
             cen_coordinates[i][2]  = coordinates[i][2] - g_z/len(coordinates)
        
        cen_coordinates = np.array(cen_coordinates)
        if self.optimal_tilt:
            #calculate the eigen values and vigen vectors
            cen_coordinates_cov=np.cov(cen_coordinates.transpose())
            #eigenvalues = np.linalg.eigvals(cen_coordinates_cov)
            #eigenvecs = np.linalg.eig(cen_coordinates_cov)
            w, v = np.linalg.eig(cen_coordinates_cov)
            finalised_cart =np.matmul(cen_coordinates,v)
        else:
            finalised_cart=deepcopy(cen_coordinates)
            
        rounded_cart=[]
        for i in range(0,len(coordinates)):
            rounded_cart.append([round(finalised_cart[i][0]),round(finalised_cart[i][1]),round(finalised_cart[i][2])])
        return rounded_cart

