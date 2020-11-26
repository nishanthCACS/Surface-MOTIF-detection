# to clear the work space
rm(list=ls())
library(bio3d)
library(plot3D)
library(plotly)
library(stringr)
library(geometry)#fort trangulation

MOTIF_initialize_func <- function(name_pdb){
# name_pdb='1m1l.pdb'
#choose the radius to MOTIF centor atom check from this radius howmany fells
boundry_radius =3.5
# how many amino acids surrounds to the center atom
neighbors=4
pdb <- read.pdb(name_pdb)
# print(pdb)
ca.inds<-atom.select(pdb,"calpha")

# plot.bio3d(pdb)
#head(pdb$xyz)
head(pdb$atom[ca.inds$atom,])
dim(pdb$atom[ca.inds$atom,])[1]
# take consider only the xyz coordinates
xyz_origin <- pdb$atom[1:dim(pdb$atom[ca.inds$atom,])[1],9:11]
# ctr+L to clear the screen

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
xyz.pca <- prcomp(xyz_origin,
                  center = TRUE,
                  scale. = FALSE) 
xyz <-predict(xyz.pca, newdata=xyz_origin)
# the xyz function change the cartesian coordinate to polar coordinate
xyz_polar <- function(xyz)
{
  # the xyz function change the cartesian coordinate to polar coordinate
  x<-xyz[1:dim(pdb$atom[ca.inds$atom,])[1],1]
  y<-xyz[1:dim(pdb$atom[ca.inds$atom,])[1],2]
  z<-xyz[1:dim(pdb$atom[ca.inds$atom,])[1],3]
  
  #since the PCA centralize
  xyz_num=cbind(x,y,z)
  
  cart2pol <- function(x, y)
  {
    r <- sqrt(x^2 + y^2)
    if (x>0){
      t <- atan(y/x)# check it is in the all region quadrant
    }
    else if(x<0 & y>0){
      t <- c(pi+atan(y/x))# check it is in the second
    }
    else{
      t <- c(-pi+atan(y/x))# check it is in the second
    }
    c(r,t)
  }
  # first create an empty matrices to hold the data
  xy_polar<-matrix(NA, nrow = dim(xyz_num)[1], ncol = 2)#here first coloumn contain the radius and the second coloumn contain angle
  for(i in 1:dim(xyz_num)[1]) {
    # # find the radius of x-y coordinate
    # xy_polar[i,1]=cart2pol(d_cal_init[i,1],d_cal_init[i,2])[1]
    # xy_polar[i,2]=cart2pol(d_cal_init[i,1],d_cal_init[i,2])[2]
    xy_polar[i,1]=cart2pol(xyz_num[i,1],xyz_num[i,2])[1]
    xy_polar[i,2]=cart2pol(xyz_num[i,1],xyz_num[i,2])[2]
  }
  
  #xy to z coordinate angle
  xy_z_polar<-matrix(NA, nrow = dim(xyz_num)[1], ncol = 2)#here first coloumn contain the radius and the second coloumn contain angle
  for(i in 1:dim(xyz_num)[1]) {
    # find the radius of xy-to z coordinate
    # #take the xy radius
    # xy_z_polar[i,1]=cart2pol(xy_polar[i,1],d_cal_init[i,3])[1]
    # xy_z_polar[i,2]=cart2pol(xy_polar[i,1],d_cal_init[i,3])[2]
    xy_z_polar[i,1]=cart2pol(xy_polar[i,1],xyz_num[i,3])[1]
    xy_z_polar[i,2]=cart2pol(xy_polar[i,1],xyz_num[i,3])[2]
  }
  
  # return(list(final_polar=cbind(xy_z_polar[1:dim(xyz)[1],1],xy_polar[1:dim(xyz)[1],2],xy_z_polar[1:dim(xyz)[1],2]),d_final=d_final))
  return(final_polar=cbind(xy_z_polar[1:dim(xyz)[1],1],xy_polar[1:dim(xyz)[1],2],xy_z_polar[1:dim(xyz)[1],2]))
  
}

final_polar<-xyz_polar(xyz)#to find out the normalized centered coordinates
# polar_function_object<-xyz_polar(xyz)#to find out the normalized centered coordinates
# final_polar= polar_function_object$final_polar
# d_final=polar_function_object$d_final

atom_angle_range <- function(angle_chk,resolution)
{
  #this function will give the range where the atom belongs
  range_xy_angle <- matrix(NA, nrow = length(angle_chk), ncol = 1)
  angle_list_rad=seq(from=-180, to=180, by=resolution)*pi/180
  for(i in 1:dim(final_polar)[1]){
    for(j in 2:length(angle_list_rad)-1){
      if(angle_list_rad[j]<angle_chk[i] & angle_chk[i]<angle_list_rad[j+1]){
        range_xy_angle[i]=j
      }
    }
  }
  range_xy_angle
}

surface_select_atom <- function(final_polar,resolution){
  # this function find out the surface atoms
  angle_xy=final_polar[1:dim(final_polar)[1],2]
  angle_xy_z=final_polar[1:dim(final_polar)[1],3]
  
  range_of_xy=atom_angle_range(angle_xy,resolution)
  range_of_xy_z=atom_angle_range(angle_xy_z,resolution)
  radius=final_polar[1:dim(final_polar)[1],1]
  # then just check the
  checked_range_element <- list();
  # selected_atom <- list();
  selected_atom<-c()
  for(i in 1:length(range_of_xy)){
    r=radius[i]
    s=i#to hold the selected item in surface
    if(is.element(i, checked_range_element)){
    }
    else{
      checked_range_element <- c(checked_range_element,i);
      for(j in 1:length(range_of_xy)){
        if(range_of_xy[i]==range_of_xy[j]){
          if(range_of_xy_z[i]==range_of_xy_z[j]){
            if(r<radius[j]){
              r=radius[j]
              s=j
            }
          }
          else{
            if(r<radius[j]){
              r=radius[j]
              s=j
            }
          }
        }
      }
      # selected_atom<- c(selected_atom,s);
      selected_atom<-rbind(selected_atom,s)
    }
  }
  unique(selected_atom)
}

resolution=1# angle degree of resolution
angle_range=seq(from=-180, to=180, by=resolution)
selected_calphas=surface_select_atom(final_polar,resolution)# this contain the surface selected atoms
m=selected_calphas
# length(m)

# the xyz function change the cartesian coordinate to polar coordinate
x<-xyz[1:dim(pdb$atom[ca.inds$atom,])[1],1]
y<-xyz[1:dim(pdb$atom[ca.inds$atom,])[1],2]
z<-xyz[1:dim(pdb$atom[ca.inds$atom,])[1],3]
# change the xyz coordinates in matrix form
xyz_num=cbind(x,y,z)
s_x=xyz_num[m,1]
s_y=xyz_num[m,2]
s_z=xyz_num[m,3]

ca.inds<-atom.select(pdb,"calpha")
dim(pdb$atom[ca.inds$atom,])[1]
# To check the residues
aa3<-pdb$atom$resid[ca.inds$atom]
aminoacids_calpha<-aa321(aa3)# change the format to one letter
# to find out the surface atoms residues -----------------------------
sur_res<-matrix(NA, nrow = length(selected_calphas), ncol = 1)
sur_res_cor<-matrix(NA, nrow = length(selected_calphas), ncol = 3)
sur_res_polar_cor<-matrix(NA, nrow = length(selected_calphas), ncol = 3)
for(i in 1:length(selected_calphas)){
  sur_res[i]=aminoacids_calpha[selected_calphas[i]]
  sur_res_cor[i,1:3]=xyz_num[selected_calphas[i],1:3]
  sur_res_polar_cor[i,1:3]=final_polar[selected_calphas[i],1:3]
  # sur_res_cor[i,1:3]= round(2*d_final[selected_calphas[i],1:3],0)
  #take shifted coordinate to [0,0,0] using i-min(i) , j-min(j), and k-min(k)
  # multiply each coordinate by 2 (just for increasing the resolution) and then round them to decimal numbers.
}

# creating results for x-y coordinates
# first create empty matrx for store the results

# then from the surface coordinates pickup the motifs
atom_distance <- function(cor_1,cor_2){
  #this function calculate the distance between the atoms
  distance=sqrt((cor_1[1]-cor_2[1])^2+(cor_1[2]-cor_2[2])^2+(cor_1[3]-cor_2[3])^2)
  return(distance)
}
# # check 
# d1=atom_distance(sur_res_cor[1,1:3],sur_res_cor[2,1:3])

#then find out the all surface atom distances
distances_sur_atoms<-matrix(NA, nrow = length(selected_calphas), ncol = length(selected_calphas))
for(i in 1:length(selected_calphas)){
  for(j in 1:length(selected_calphas)){
    distances_sur_atoms[i,j]=atom_distance(sur_res_cor[i,1:3],sur_res_cor[j,1:3])
  }
}
return(list(sur_res_cor = sur_res_cor, sur_res_polar_cor = sur_res_polar_cor, sur_res = sur_res))
}
 
getwd()
# list.files()
pdb_ids_obj <- read.table("unique_tsg.csv",sep=",", header=FALSE)
pdb_ids_obj_2<-pdb_ids_obj[1:dim(pdb_ids_obj)[1],1]
pdb_ids<-matrix(NA, nrow = length(pdb_ids_obj_2), ncol = 1, byrow = TRUE)
for (i in 1:length(pdb_ids_obj_2)){
  pdb_ids[i]<-toString(pdb_ids_obj_2[i])
}
pdb_ids_naming<-matrix(NA, nrow = length(pdb_ids_obj_2), ncol = 1, byrow = TRUE)
for (i in 1:length(pdb_ids_obj_2)){
  pdb_ids_naming[i]<-toString(pdb_ids_obj_2[i])
}

for(m in  13:length(pdb_ids_naming)){
  print(paste("Id number in progress: ",m))
  #to find out the normalized centered coordinates
  MOTIF_INI_function_object<-MOTIF_initialize_func(pdb_ids_naming[m])
  sur_res_cor = MOTIF_INI_function_object$sur_res_cor
  sur_res_polar_cor = MOTIF_INI_function_object$sur_res_polar_cor
  sur_res = MOTIF_INI_function_object$sur_res
  
  write.table(sur_res_cor, file=str_c(pdb_ids_naming[m],"_surf_atoms.csv"),sep = ",",row.names = F)
  write.table(sur_res_polar_cor, file=str_c(pdb_ids_naming[m],"_surf_polar_atoms.csv"),sep = ",",row.names = F)
  write.table(sur_res, file=str_c(pdb_ids_naming[m],"_aminoacids.csv"),sep = ",",row.names = F)
  
  
  print(paste('------------------------done: ',pdb_ids_naming[m]))
}