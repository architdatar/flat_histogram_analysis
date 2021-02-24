#%%
"""
Contains methods for analyzing the isotherms and drawing insights. 
"""

import scipy.interpolate as sci
import scipy.optimize as sco
from functools import reduce 

from . import plt, np, pd, scipy

class IsothermAnalysis():
    """
    Analyze isotherms and draw insights from them. 

    Methods:
    --------
    convert_to_molecules_per_framework: Isotherm in different units convert to molecules 
        per framework. Currently, only "kg/kg" is supported. 
    convert_from_molecules_per_framework: Isotherm in molecules per framework 
        different units convert to other units. Currently, only "kg/kg" is supported.
    get_step_pressure: Get pressure at isotherm step (at loading at a 
        fraction of geometric loading).
    """
    def __init__(self):
        """
        """
        self.WATER_MASS = 18 #g/mol

    def convert_to_molecules_per_framework(self, frame_mass, input_unit="kg/kg"):
        """
        Returns conversion factor for units to molecules / framework. 

        Other units will be added later. 
        Parameters
        ------------
        frame_mass: Mass of framework in g/mol.
        Returns
        ------------
        Conversion factor to convert from current unit to molecule/framework
        """
        if input_unit == "kg/kg":
            return frame_mass / self.WATER_MASS
    
    def convert_from_molecules_per_framework(self, frame_mass=None, 
                                            output_unit="kg/kg"):
        """
        Returns conversion factor for units from molecules / framework
        to desired units. 

        Other units will be added later. 
        Parameters
        ------------
        frame_mass (optional): Mass of framework in g/mol.
        Returns
        ------------
        Conversion factor to convert from molecule/framework to current unit.
        """
        if output_unit == "kg/kg":
            return self.WATER_MASS / frame_mass

    def get_step_pressure(self, pressure, loading, step_loading):
        """
        Get pressure of step from the isotherm. 

        Parameters
        ------------
        pressure: Isotherm pressure values.
        loading: Isotherm loading values.
        step_loading: Loading corresponding to the step as a fraction of geometric loading. 

        Returns
        ------------
        step_pressure: Pressure at step of the isotherm.
        """
        interp_y = sci.interp1d(loading, pressure) 
        try:
            step_p =  interp_y(step_loading)
            step_p = np.float(step_p)
            #print(step_p)        
            if np.isnan(step_p):
                step_pressure = 10000
            else:
                step_pressure = step_p
        except:
            step_pressure=10000
        return step_pressure


if __name__ == "__main__":
    iso=IsothermAnalysis()
    df_iso = pd.read_csv("NAVFUM_fullCoRE_298_918_spline_C-iso-5pt_reweighted.txt", 
                        header=0)
    step_pressure = iso.get_step_pressure(df_iso["Pressure"], df_iso["NLoading"], 150)
    print(step_pressure)


# %%
