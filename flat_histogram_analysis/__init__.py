"""Importing the module files and related libraries. 
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
