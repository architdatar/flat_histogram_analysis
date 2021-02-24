"""
Scan through the C data available. Make isotherms, get location of step. 
"""
#%%
#%load_ext autoreload
%reload_ext autoreload
%autoreload 2

"""
A use case for the package is illustrated. The goal
is to analyze the C data for a structure and generate
an adsorption isotherm. Then, using that isotherm, 
get the pressure at the step. 
"""

import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import time

#import flat_histogram_analysis
from flat_histogram_analysis.NVTW_data_analysis import NVTW_analysis
from flat_histogram_analysis.isotherm_analysis import IsothermAnalysis

wh = NVTW_analysis()
iso = IsothermAnalysis()

water_mass=iso.WATER_MASS #provided in g/mol

P0_dict={'TIP4P': {'298': 4594, '338': 35712}} #P_sat values. 
P_ref_dict = {'TIP4P': {'298': 918, '338': 5357}} #Reference pressure at which 
            #simulations were actually done. Certain percent of Psat. 

filename_root = "C-iso-5pt" 
water_model = "TIP4P"

#attribute_dataframe.
df_summary = pd.DataFrame()

time_start_og = time.time()

hydrophobic_str = pd.read_csv("~/Downloads/hydrophobic_str", 
            usecols=[1], names=["Name_full"], skiprows=1)

name_root = "NAVFUM" #NAVFUM is the CSD refcode of the MOF structure.
name = f"{name_root}_fullCoRE" 

frame_mass = 41433 #g/mol
#geo_load =  7.9 #molecules / framework. 
#unit_cells = np.array([4,4,2]) #Based on the CIF file to ensure periodicity.
#geo_load_kgkg = geo_load * water_mass / frame_mass * unit_cells.prod() #geometric water loading in kg/kg
geo_load_kgkg =  0.1098 #kg/kg # can be obatined by multiplying pore volume * density of water. 

print(f"Trying for {name}")
#if True:
try:
    df_original = pd.read_table(f"{name}_fCoRE.dat", 
                                delim_whitespace=True, header=0)

    # Output paths
    C_image_filepath =  f"water_model_C_image_{name}_" \
                                    f"{filename_root}.png"
    P_image_filepath = f"water_model_P_image_{name}_" \
                                    f"{filename_root}.png"
    energy_image_filepath = f"water_model_energy_image_{name}_" \
                                    f"{filename_root}.png"
    isotherm_figure_filepath = f"water_model_isotherm_image_{name}_" \
                                            f"{filename_root}.png"

    #Initialize the figure objects. 
    fig_C = wh.plot_C_values(fig=None, figure_title=name,
                            plot_data=False, final_adjustment=False)
    fig_P = wh.plot_P_values(fig=None, figure_title=name,
                            plot_data=False, final_adjustment=False)
    fig_energy = wh.plot_energy_values(fig=None, figure_title=name,
                            plot_data=False, final_adjustment=False)

    #Initializing isotherm plots. The idea was to keep them more flexible
    #then the other plots, but in later versions, it could be written as a method. 
    [fig_iso, [ax_mol, ax_kgkg] ]= plt.subplots(nrows=1, ncols=2, figsize=(14,6))
    ax_mol.set_xlabel("Fugacity (Pa)")
    ax_mol.set_ylabel("Loading (mol/framework)")
    ax_mol.set_title(name)
    ax_kgkg.set_xlabel("PRel")
    ax_kgkg.set_ylabel("Loading (kg/kg)")
    ax_kgkg.set_title("Isotherms")
    
    step_pressure_dict = {} #Inititalize this for every structure. 

    for T_ref_index, T_ref in enumerate([298, 338]):
        P_sat = P0_dict[water_model][f"{T_ref}"]
        P_ref = P_ref_dict[water_model][f"{T_ref}"]
        
        pressure_list = np.concatenate([np.linspace(0, P_sat, 30, endpoint=False), 
                        np.linspace(P_sat, 20*P_sat, 20)])

        df = df_original[df_original["Temp"]==T_ref].copy(deep=True)

        #interp_method="regression"
        interp_method = "spline"

        #This is to calculate isotherms where where pressure list extends beyond Prel=1.0.
        C_data_filepath = f"{name}_{T_ref:d}_{P_ref:.0f}_" \
                                        f"{interp_method}_{filename_root}_C-Data.xlsx"
        full_isotherm_filepath = f"{name}_{T_ref}_{P_ref:.0f}_" \
                                        f"{interp_method}_{filename_root}_reweighted.txt" #Name when we are 
                                                                                #interested in going beyond P_rel =1. 

        #Calculates C data.
        [df, df_C] = wh.interpolate_C_data( df, name, T_ref, P_ref, interp_method=interp_method, filename_root=filename_root)
        [df, df_C] = wh.calculate_P_values(df, df_C, save_C_tables=False)
                
        #Plot these. 
        fig_C = wh.plot_C_values(df, df_C, fig=fig_C, plot_data=True, 
                                plot_color=plt.cm.Set1(T_ref_index), 
                                plot_label=f"T={T_ref}", final_adjustment=False,
                                save_image=True, image_filepath=C_image_filepath)
        fig_P = wh.plot_P_values(df, df_C, fig=fig_P, plot_data=True, 
                                plot_color=plt.cm.Set1(T_ref_index), 
                                plot_label=f"T={T_ref}", final_adjustment=False, 
                                save_image=True, image_filepath=P_image_filepath)
        fig_energy = wh.plot_energy_values(df, df_C, fig=fig_energy, plot_data=True,
                                plot_color=plt.cm.Set1(T_ref_index), 
                                plot_label=f"T={T_ref}", final_adjustment=False, 
                                save_image=True, image_filepath=energy_image_filepath)
        
        #Get isotherm data. 
        [df_iso, df_C_iso] = wh.get_isotherm_data(df_C, temp_list=[T_ref], 
                                    pressure_list=pressure_list, 
                                    P_ref=P_ref, T_ref=T_ref, 
                                    interp_method=interp_method, 
                                    filename_root=filename_root, 
                                    save_isotherm_file=False,
                                    #isotherm_filepath=isotherm_filepath, #This file
                                        # only stores pressure and loading. But, we also want 
                                        # relative pressure and loading in kg/kg, so we save another file. 
                                    save_full_C_data_file=True,
                                    C_data_filepath=C_data_filepath)

        df_iso.loc[df_iso["Pressure"]==0, "NLoading"] = 0 #Set pressure 0 to loading 0.

        df_iso["PRel"] = df_iso["Pressure"] / P_sat
        df_iso["Loading_kgkg"] = df_iso["NLoading"] * \
                        iso.convert_from_molecules_per_framework(frame_mass, output_unit="kg/kg")
        df_iso.to_csv(full_isotherm_filepath, index=False)

        #Add to isotherm plots. 
        ax_mol = wh.plot_isotherm_figure(df_iso, ax=ax_mol, plot_data=True,
                        x_value="Pressure", y_value="NLoading", 
                        plot_color=plt.cm.Set1(T_ref_index),
                        plot_label=f"T={T_ref}", final_adjustment=False)
        ax_kgkg = wh.plot_isotherm_figure(df_iso, ax=ax_kgkg, plot_data=True,
                        x_value="PRel", y_value="Loading_kgkg", 
                        plot_color=plt.cm.Set1(T_ref_index),
                        plot_label=f"T={T_ref}", final_adjustment=True)

        step_pressure = iso.get_step_pressure(df_iso["Pressure"], df_iso["Loading_kgkg"], geo_load_kgkg/2)
        step_p_rel = iso.get_step_pressure(df_iso["PRel"], df_iso["Loading_kgkg"], geo_load_kgkg/2)

        #Save step pressure data into a dictionary. 
        step_pressure_dict[f"{T_ref}"] = {"step_pressure": step_pressure, "step_prel": step_p_rel, "step_loading": geo_load_kgkg/2} 

    #Plot the geometric loading lines and annotate the step pressures.
    # For T=298 and 338 K referenced explicitly.  
    geo_load_mol = geo_load * unit_cells.prod()
    ax_mol.set_xlim(ax_mol.get_xlim())
    ax_mol.set_ylim(ax_mol.get_ylim())
    ax_mol.hlines(geo_load_mol, ax_mol.get_xlim()[0], ax_mol.get_xlim()[1])
    ax_mol.hlines(geo_load_mol / 2, ax_mol.get_xlim()[0], step_pressure_dict["298"]["step_pressure"], linestyles='--')
    ax_mol.vlines(step_pressure_dict["298"]["step_pressure"], ax_mol.get_ylim()[0], geo_load_mol/2, linestyles='--')
    ax_mol.hlines(geo_load_mol / 2, ax_mol.get_xlim()[0], step_pressure_dict["338"]["step_pressure"], linestyles='--')
    ax_mol.vlines(step_pressure_dict["338"]["step_pressure"], ax_mol.get_ylim()[0], geo_load_mol/2, linestyles='--')

    ax_kgkg.set_xlim(ax_kgkg.get_xlim()) #this is for the expanded structures. 
    ax_kgkg.set_ylim(ax_kgkg.get_ylim())
    ax_kgkg.hlines(geo_load_kgkg, ax_kgkg.get_xlim()[0], ax_kgkg.get_xlim()[1])
    ax_kgkg.hlines(geo_load_kgkg / 2, ax_kgkg.get_xlim()[0], step_pressure_dict["298"]["step_prel"], linestyles='--')
    ax_kgkg.vlines(step_pressure_dict["298"]["step_prel"], ax_kgkg.get_ylim()[0], geo_load_kgkg/2, linestyles='--')
    ax_kgkg.hlines(geo_load_kgkg / 2, ax_kgkg.get_xlim()[0], step_pressure_dict["338"]["step_prel"], linestyles='--')
    ax_kgkg.vlines(step_pressure_dict["338"]["step_prel"], ax_kgkg.get_ylim()[0], geo_load_kgkg/2, linestyles='--')
        
    #Close the figures. 
    wh.plot_C_values(fig=fig_C, 
                            plot_data=False, final_adjustment=True, 
                            save_image=True, image_filepath=C_image_filepath,
                            close_fig=True)
    wh.plot_P_values(fig=fig_P, 
                            plot_data=False, final_adjustment=True,
                            save_image=True, image_filepath=P_image_filepath, 
                            close_fig=True)
    wh.plot_energy_values(fig=fig_energy, 
                            plot_data=False, final_adjustment=True, 
                            save_image=True, image_filepath=energy_image_filepath,
                            close_fig=True)
    #Isotherm figure can be closed separately. 
    fig_iso.savefig(isotherm_figure_filepath, bbox_inches="tight", dpi=300)
    plt.close("all")

    #Get all the relevant information out for the structure in this dictionary.
    property_dict={"Name": name,  "Geo_load": geo_load_kgkg, 
            "StepPressure_298": step_pressure_dict["298"]["step_pressure"],
            "StepPRel_298": step_pressure_dict["298"]["step_prel"],
            "StepPressure_338": step_pressure_dict["338"]["step_pressure"],
            "StepPRel_338": step_pressure_dict["338"]["step_prel"]  }
    df_summary = df_summary.append(property_dict, ignore_index=True)
        
except Exception as e:
#else:
    #e=""
    print (f"Problem {name} {e}")
    #continue

cols = df_summary.columns.tolist()
cols.remove("Name") #Rearranging columns to ensure that name is first.
cols.insert(0, "Name")
df_summary = df_summary[cols]
display(df_summary)

time_end_og = time.time()
time_taken_og =  time_end_og - time_start_og
print(f"Time taken overall: {time_taken_og}")

df_summary.to_excel("StepPressure.xlsx", 
                   index=False) 

# %%
