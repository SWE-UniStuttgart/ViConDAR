# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 15:08:09 2019

@author: pettas
"""
import matplotlib.pyplot as plt  # matplotlib for some plotting
import numpy as np  # numeric python functions
import pandas as pd  # need this to load our data from the csv files

from pyconturb import gen_turb, gen_spat_grid  # generate turbulence, useful helper
from pyconturb.sig_models import iec_sig  # IEC 61400-1 turbulence std dev
from pyconturb.spectral_models import kaimal_spectrum  # Kaimal spectrum
from pyconturb.wind_profiles import constant_profile, power_profile  # wind-speed profile functions
from pyconturb.wind_profiles import get_wsp_values,constant_profile, power_profile, data_profile
from pyconturb.spectral_models import kaimal_spectrum, data_spectrum
from pyconturb.sig_models import iec_sig, data_sig # IEC 61400-1 turbulence std dev



y = np.linspace(-90, 90, 41)  # 11 lateral points from -50 to 50 (center @ 0)
z = np.linspace(30, 210, 41)  # 13 vertical points from 40 to 160 (center @ 100)
spat_df = gen_spat_grid(y, z)  # if `comps` not passed in, assumes all 3 components are wanted
wsp_func = power_profile
sig_func = iec_sig


kwargs = {'u_ref': 10, 'z_ref': 120, 'alpha':0.3,'T': 600, 'dt': 2}
turb_df = gen_turb(spat_df,wsp_func=wsp_func, **kwargs)

turb_df.to_csv('DTU10MW_V10_SH03')