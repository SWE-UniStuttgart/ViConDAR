# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 13:49:11 2020

@author: Davide Conti davcon@dtu.dk

This script serves to read binary Mann turbulence boxes (u,v,w) and to 
reshape the boxes in order to be compliant with the required input fields to
the Vicondar framework. 
"""

import numpy as np
import os
import scipy
from scipy import io
from scipy.io import savemat


    
############################################################################
###### Load a turbulence box generated with the Mann.exe 
############################################################################

filename = r'./turb1001_%s.bin' #defines the name convention for the u,v,w files from the mann generator
output_name = './vasilis.mat' # defines the output name of the .mat file for vicondar



# turbulence box dimensions (same inputs used in the Mann.exe)
nx,ny,nz = 8192, 32,32 
dy, dz = 6.5,6.5, 
dx= 0.02
ly, lz = dy*ny,dz*nz, 
Tsim = 600 # simulation length in seconds. 

# Load the binary files of the turbulence boxes by the Mann.exe (these are residual fields)
u_mann, v_mann, w_mann = [np.fromfile(filename % uvw, np.dtype('<f'), -1).reshape(nx , ny, nz) for uvw in ['u', 'v', 'w']]        

# Ambient wind speed (m/s), a vertical shear profile can be used 
URef = 10 

# Spatial grid of the turbulence box
Y_coord= np.arange(-ly/2, ly/2, ly/ny)
Z_coord= np.arange(-lz/2, lz/2, lz/nz)

# Sampling frequency/period of the wind field in turbulence box. 
fs = nx/Tsim
dt = 1/fs
    

###########################################################################
###  Convert the wind fields (u,v,w) to the Vicondar's convetion 
###########################################################################

# The symmetry of the wind field is required in Vicondar (that's why we have +1)
u_vicondar = np.zeros(shape=(nx,ny+1,nz+1)) 
v_vicondar = np.zeros_like(u_vicondar)
w_vicondar = np.zeros_like(u_vicondar)
for i in range(ny+1):
    for j in range(nz+1):
        if (i == ny) & (j < nz) :
            u_vicondar[:,i,j] = u_mann[:,i-1,j]
            v_vicondar[:,i,j] = v_mann[:,i-1,j]
            w_vicondar[:,i,j] = w_mann[:,i-1,j]
        elif (j == nz) & (i <ny) :
            u_vicondar[:,i,j] = u_mann[:,i,j-1]
            v_vicondar[:,i,j] = v_mann[:,i,j-1]
            w_vicondar[:,i,j] = w_mann[:,i,j-1]
        elif (i == ny) & (j == nz):
            u_vicondar[:,i,j] = u_mann[:,i-1,j-1]                
            v_vicondar[:,i,j] = v_mann[:,i-1,j-1]                
            w_vicondar[:,i,j] = w_mann[:,i-1,j-1]                
        else:             
            u_vicondar[:,i,j] = u_mann[:,i,j]
            v_vicondar[:,i,j] = v_mann[:,i,j]
            w_vicondar[:,i,j] = w_mann[:,i,j]
 
#### Derive parameters compliant with Vicondar (i.e. name convention, coordinate systems)
   
# Required in Vicondar
T_offset = 6.125  
# Box coordinates
nt = nx; 
y_hori = np.hstack((Y_coord,-Y_coord[0]))
z_vert = np.hstack((Z_coord,-Z_coord[0]))
t, y1, z1 = np.arange(0,nx*dt,dt).astype('d'), y_hori.astype('d'), z_vert.astype('d')
Y1,Z1 = np.meshgrid(y1,z1)

# Reshape u,v,w arrays 
u_VCD =np.zeros(shape=(int(ny)+1,int(nt),int(nz)+1))
v_VCD =np.zeros_like(u_VCD); w_VCD =np.zeros_like(u_VCD)
for i in range(nz+1):
    u_VCD[i,:,:] = u_vicondar[:,:,i] 
    v_VCD[i,:,:] = v_vicondar[:,:,i] 
    w_VCD[i,:,:] = w_vicondar[:,:,i] 
    
# Convert data as double type (compatibility with Matlab)
u_VCD,v_VCD,w_VCD = u_VCD.astype('d'), v_VCD.astype('d'), w_VCD.astype('d') 
Y1,Z1 = Y1.astype('d'),Z1.astype('d')

# Save outputs in windfield structure
nt = float(nt);  ny = float(ny)+1; URef = float(URef)
nz = float(nz) +1;
# Create the "windfield" structure to Vicondar 
grid = {'nt': nt, 'ny': ny, 'nz': nz, 'dt': dt, 'dy': dy, 'dz': dz,'t': t, 'y': y1, 'z': z1,'Y': Y1, 'Z': Z1}
windfield = {'grid': grid, 'ny': ny, 'nz': nz, 'dt': dt,'u': u_VCD,'v': v_VCD,'w': w_VCD, 'URef': URef,'T_offset':T_offset}
windfield= {'windfield': windfield}

# Create a wind field variable     
ID_file = 1 # random id. 
scipy.io.savemat(r'%s' %(output_name), windfield)    

    




