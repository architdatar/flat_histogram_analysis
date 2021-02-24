"""
flat_histogram_analysis 
========================

Analyzes the data from the flat histogram Monte Carlo simulations 
for adsorption using the NVT+W method. Largely used to develop the calculations
used in the JPCC paper (DOI: 10.1021/acs.jpcc.0c11082). 

Basically, the data from the simulations--the C data--is used to get
the relative weights of each loading macrostate and then used to 
generate adsorption isotherms which are then analyzed. 

Modules:
--------
NVTW_data_analysis: Contains methods for interpolating the C data, 
    creating relative weights, and reweighting to generate isotherm. 
isotherm_analysis: Analyzes the generated adsorption isotherms and 
    calculates their interesting properties such as adsorption 
    isotherm step. 
"""
#%%
import os
import sys
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt

#from .NVTW_data_analysis import NVTW_analysis
#from .isotherm_analysis import IsothermAnalysis
import time

pd.set_option('display.max_rows', 800)
pd.set_option('display.expand_frame_repr', False)

#print(f"CWD: {os.getcwd()}")
dir_path = os.path.dirname(os.path.realpath(__file__))
#print(f"File path: {dir_path}")
plt.style.use(os.path.join(dir_path, 'mypaper_arial.mplstyle'))


# %%
