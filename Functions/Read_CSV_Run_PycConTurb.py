# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 19:02:30 2019

@author: fcosta
"""
"""
Created on Wed Jul 10 09:30:19 2019

@author: fcosta

Read_CSV_Run_PyConturb:
    
Script that imports outputs from ViConDAR and creates a turbulence Wind field using PyConTurb. It is meant to be called directly 
from matlab with arguments defined in the matlab code in order to ryn PyConTurb from Matlab. Expects that python is in the root and 
the base environment is compatible with pyconturb.
The output of this script is a turbulence wind field saved in .csv format.
I/O:
    The main inputs are choosen in matlab.
    - Variables contains info about:
        > Total time of the time series [s]
        > time step of the time series [s]
        > reference mean velocity [m/s]
        > zhub [m]
        > coherence, wsp, sig and spec models
        > seed
        > interpolation data
        > Turbulence class
    - GridY and gridZ contains info about the used grid of points in space (velocities[m/s])
    - con_df contains info about the velocities of the constrained points along "Total time" in steps of "dt" seconds
    - PyConTurb_Simulation: Output for ViConDAR in .csv format.

"""
#%% Import

import sys 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  # matplotlib for some plotting
import os
from pyconturb.wind_profiles import get_wsp_values,constant_profile, power_profile, data_profile
from pyconturb.spectral_models import kaimal_spectrum, data_spectrum
from pyconturb.sig_models import iec_sig, data_sig # IEC 61400-1 turbulence std dev
from pyconturb import gen_turb, gen_spat_grid, TimeConstraint # generate turbulence, useful helper

# parse arguments
condf_name = sys.argv[1]
gridy_name = sys.argv[2]
gridz_name = sys.argv[3]
variables_name = sys.argv[4]
savename   = sys.argv[5]

#  Reading CSV File from Matlab an running PyConTurb

Variables = pd.read_csv(variables_name,index_col=None)
GridY = pd.read_csv(gridy_name,index_col=None)
GridZ = pd.read_csv(gridz_name,index_col=None)
#import pdb; pdb.set_trace()
con_df=pd.read_csv(condf_name,index_col=0)

#import pdb; pdb.set_trace()

## Time Constraints:
con_tc =  TimeConstraint(con_df)  # Constraining in time
con_tc.index = con_tc.index.map(lambda x: float(x) if (x not in 'kxyz') else x)  
                 
    # Inputs to constrained Turbulence:
spat_df = gen_spat_grid(GridY.iloc[0, :], GridZ.iloc[0, :])  # create our spatial pandas dataframe with our grid data.
kwargs  = {'u_ref': Variables.iloc[0, 2], 'z_ref': Variables.iloc[0, 4],'z_hub': Variables.iloc[0, 4], # necessary keyword arguments for IEC turbulence
              'T': con_tc.get_T(), 'dt': con_tc.get_time().index[1],'turb_class': Variables.iloc[0,-1],'alpha':Variables.iloc[0,-2],'seed':Variables.iloc[0,9]}  # simulation length (s) and time step (s)

wsp_funcin   = Variables.iloc[0,6]
sig_funcin   = Variables.iloc[0,7]
spec_funcin  = Variables.iloc[0,8]
interp_data2 = Variables.iloc[0,10]
nf_chunkin   = Variables.iloc[0,3]     
coh_model_in = Variables.iloc[0, 5]  

if interp_data2=="Take_list":
    # here we have to define all the cases in wsp_func, spec_func and sig_func that we want to use as inputs for "gen_turb":
    if wsp_funcin=='power_profile':
        wsp_func = power_profile
    elif wsp_funcin=='constant_profile':
        wsp_func = constant_profile
    elif wsp_funcin=='data_profile':
        wsp_func = data_profile
#        wsp_func(GridY.iloc[0, :],GridZ.iloc[0, :],con_tc=None,**kwargs)                  
    if sig_funcin=='iec_sig':
       sig_func = iec_sig  
    elif sig_funcin=='data_sig':
       sig_func = data_sig       
#       sig_func(0,GridY.iloc[0, :],GridZ.iloc[0, :],**kwargs)      data_spectrum should be used            
    if spec_funcin=='kaimal_spectrum':
       spec_func = kaimal_spectrum
    elif spec_funcin=='data_spectrum':
       spec_func = data_spectrum        
    interp_data=[wsp_func,sig_func,spec_func]                   
elif interp_data2=='none' or interp_data2=='None':
     interp_data = 'none'
elif interp_data2=='all':
     interp_data = 'all'
 
#    import pdb; pdb.set_trace()            
# Executing PyConTurb
if interp_data2=="Take_list":
    # import pdb; pdb.set_trace()
#        sim_turb_df = gen_turb(spat_df,con_tc=con_tc, wsp_func=interp_data[0],sig_func=interp_data[1],spec_func=interp_data[2],coh_model_in=Variables.iloc[0, 5],**kwargs)
    sim_turb_df = gen_turb(spat_df,con_tc=con_tc, wsp_func=interp_data[0],sig_func=interp_data[1],spec_func=interp_data[2],nf_chunk= nf_chunkin,verbose=True,**kwargs)
else:
    sim_turb_df = gen_turb(spat_df,con_tc=con_tc,interp_data=interp_data, **kwargs)

#  Save the turbulent Windfield in .csv format;
sim_turb_df.to_csv(savename)