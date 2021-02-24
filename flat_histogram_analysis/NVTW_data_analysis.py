"""
Contains methods for interpolating the C data, 
creating relative weights, and reweighting to generate isotherm. 
"""
import scipy.interpolate as sci
import scipy.optimize as sco
from functools import reduce 

#Defined in __init__.py. The plt class has been modified there.
from . import plt, np, pd, scipy

class NVTW_analysis():
    """
    Performs the analysis on the data from the NVT+W method results. 

    Methods
    -------
    interpolate_C_data: Interpolates the C data to create values at every N.
    calculate_P_values: Computes the P values, along with the 
        $\Pi$ (weight) values.
    reweight_iso: Performs histogram reweighting to get the weight distribution at 
        any T, P. 
    get_isotherm_data: Uses the weight distribution to get the weighted N value at 
        a given set of T, P conditions to generate an isotherm. 
    plot_C_values: Makes a plot of the C values.
    plot_P_values: Makes a plot of the P values.
    plot_energy_values: Makes a plot of the energy moments collected from the simulation. 
    plot_isotherm_figure: Plots the isotherm figure. 
    """
    def __init__(self):
        """
        """

    def interpolate_C_data(self, df, name, T_ref, P_ref, interp_method="regression", filename_root="rd2",
                           save_C_tables=False, filepath_og="tmp_C_table_original.xlsx", filepath_C="tmp_C_table_interpolated.xlsx"):
        """
        Interpolate C data. 
        
        Parameters
        ------------
        df: Data frame containing the C data in the columns C_N-N-1, C_N-N and C_N-N+1. 
        name: Name of the structure for which the analysis is being done. 
        T_ref: Temperature at which the simulation was done. 
        P_ref: Pressure (fugacity) at which the simulation was done. 
        interp_method: Interpolation method for the C data. Either cubic "spline" or
            polynomial "regression" (default: "regression").
        filename_root = This is to identify the files resulting from the analysis (default: rd2).
        save_C_tables: Whether to save the output interpolated C data (default: False)
        filepath_og: File path to which the original C data are stored (default: tmp_C_table_original.xlsx)  
        filepath_C: File path to which the interpolated C data are stored (default: tmp_C_table_interpolated.xlsx)

        Returns
        -------
        df: The original data frame. 
        df_C: The dataframe containing the interpolated data at each N. 
        """
        
        if df[np.isnan(df['C_N-N-1'])].shape[0] > 0:
            print ("There are NaN values in the dataframe")
            print (df[np.isnan(df['C_N-N-1'])])
        
        df=df[~np.isnan(df['C_N-N-1'])] 
        
        df = df[df['N']>0]
        #df = df[df['N']>=0]
            
        #Here, we try out the ideas of using an average of C values to interpolate between the various C values, to ensure lower noise.
        df['C_N-N-1_og'] = df['C_N-N-1'] 
        df['C_N-N_og'] = df['C_N-N'] 
        df['C_N-N+1_og'] = df['C_N-N+1'] 
        df['C_N-N-1'] = df['C_N-N-1'] / df['Samples']
        df['C_N-N'] = df['C_N-N'] / df['Samples']
        df['C_N-N+1'] = df['C_N-N+1'] / df['Samples']
        
        df['Energy_mean'] = df['Energy'] / df['Samples']
        df['Energy2_mean'] = df['Energy2'] / df['Samples']
        df['Energy3_mean'] = df['Energy3'] / df['Samples']

        spacing=1 #Keep this to 1, strictly!
        
        points_1 = np.arange(df['N'].values[0], df['N'].values[-1] + 1, spacing)
        points_2 = np.arange(df['N'].values[0], df['N'].values[-1] + 1, spacing)
        points_3 = np.arange(df['N'].values[0], df['N'].values[-1] + 1, spacing)
        #points_3 = np.arange(df['N'].values[0], 220 +1 , spacing) #449 #492
        
        points_4 = np.arange(df['N'].values[0], df['N'].values[-1] + 1, spacing)
        points_5 = np.arange(df['N'].values[0], df['N'].values[-1] + 1, spacing)
        points_6 = np.arange(0, df['N'].values[-1] + 1, spacing)
        
        #df_interp_1 = df[df.index.isin(points_1-1)]
        #df_interp_2 = df[df.index.isin(points_2-1)]
        #df_interp_3 = df[df.index.isin(points_3-1)]
        
        df_interp_1 = df[df['N'].isin(points_1)]
        #df_interp_1 = df[df['C_N-N-1'] > 0.1 ]
        df_interp_2 = df[df['N'].isin(points_2)]
        df_interp_3 = df[df['N'].isin(points_3)]
        
        df_interp_4 = df[df['N'].isin(points_4)]
        df_interp_5 = df[df['N'].isin(points_5)]
        df_interp_6 = df[df['N'].isin(points_6)]
        
        #x_range=np.linspace(0, 280, 100)
        x_range=np.arange(0, df['N'].values[-1])
        
        if interp_method=="spline":
            #cubic splines:
            """
            sp1 = sci.CubicSpline(df_interp['N'], df_interp['C_N-N-1'])
            sp2 = sci.CubicSpline(df_interp['N'], df_interp['C_N-N'])
            sp3 = sci.CubicSpline(df_interp['N'], df_interp['C_N-N+1'])
            """
    
            #interpolation with log scale. 
            sp1 = sci.CubicSpline(df_interp_1['N'], scipy.log10(df_interp_1['C_N-N-1']))
            sp2 = sci.CubicSpline(df_interp_2['N'], scipy.log10(df_interp_2['C_N-N']))
            sp3 = sci.CubicSpline(df_interp_3['N'], scipy.log10(df_interp_3['C_N-N+1']))
            
            sp4 = sci.CubicSpline(df_interp_1['N'], df_interp_1['Energy_mean'])
            sp5 = sci.CubicSpline(df_interp_2['N'], df_interp_2['Energy2_mean'] )
            sp6 = sci.CubicSpline(df_interp_3['N'], df_interp_3['Energy3_mean'] )
            
            s1 = sp1(x_range)
            s2 = sp2(x_range)
            s3 = sp3(x_range)
            
            s1 = 10**(s1)
            s2 = 10**(s2)
            s3 = 10**(s3)
            
            s4 = sp4(x_range)
            s5 = sp5(x_range)
            s6 = sp6(x_range)
            
        elif interp_method=="regression":
            #Regression with polynomials.
            def polynomial(x, params): 
                cumulative=0
                for index, value in enumerate(params):
                    cumulative += value * x**index 
                return cumulative
            
            def error(params, x, y): 
                return  (polynomial(x, params) - y)
            
            def log_polynomial(x, params):
                return scipy.log10(polynomial(x, params))
            
            def log_error(params, x, y): 
                return  (log_polynomial(x, params) - y)
            
            #params_init_1 = [1,1,1,1]
            #params_init_2 = [1,1,1]
            #params_init_3 = [1,1,1]
            
            params_init_1 = [1,1,1,1]
            params_init_2 = [1,1,1]
            params_init_3 = [1,1,1,1]
            
            params_init_4 = [1,1,1]
            params_init_5 = [1,1,1]
            params_init_6 = [1,1,1]
            
            #res_1 = sco.leastsq(error, params_init_1, args=(df_interp_1['N'], df_interp_1['C_N-N-1']) , full_output=1)
            #res_2 = sco.leastsq(error, params_init_2, args=(df_interp_2['N'], df_interp_2['C_N-N']) , full_output=1)
            #res_3 = sco.leastsq(error, params_init_3, args=(df_interp_3['N'], df_interp_3['C_N-N+1']) , full_output=1)
            
            res_1 = sco.leastsq(error, params_init_1, args=(df_interp_1['N'], scipy.log10(df_interp_1['C_N-N-1'])) , full_output=1)
            res_2 = sco.leastsq(error, params_init_2, args=(df_interp_2['N'], scipy.log10(df_interp_2['C_N-N'])) , full_output=1)
            res_3 = sco.leastsq(error, params_init_3, args=(df_interp_3['N'], scipy.log10(df_interp_3['C_N-N+1'])) , full_output=1)
            
            res_4 = sco.leastsq(error, params_init_4, args=(df_interp_1['N'], df_interp_1['Energy_mean']) , full_output=1)
            res_5 = sco.leastsq(error, params_init_5, args=(df_interp_2['N'], df_interp_2['Energy2_mean']) , full_output=1)
            #res_6 = sco.leastsq(error, params_init_6, args=(df_interp_6['N'], df_interp_6['Energy3_mean']) , full_output=1)
            res_6 = sco.leastsq(error, params_init_6, args=(df_interp_6['N'], scipy.log10(-df_interp_6['Energy3_mean'])) , full_output=1)
            #res_6 = sco.leastsq(error, params_init_6, args=(df_interp_3['N'], df_interp_3['Energy3_mean']) , full_output=1)
            
            params_1=res_1[0]
            params_2=res_2[0]
            params_3=res_3[0]
            
            params_4=res_4[0]
            params_5=res_5[0]
            params_6=res_6[0]
            
            s1=polynomial(x_range, params_1)
            s2=polynomial(x_range, params_2)
            s3=polynomial(x_range, params_3)
            
            s4=polynomial(x_range, params_4)
            s5=polynomial(x_range, params_5)
            s6=polynomial(x_range, params_6)
             
            s1 = 10**(s1)
            s2 = 10**(s2)
            s3 = 10**(s3)
            
            #s4 = s4
            #s5 = s5
            s6 = -10**(s6)
                
        df_C = pd.DataFrame(columns=['N', 'C_N-N-1', 'C_N-N', 'C_N-N+1'])
        df_C['N'] = x_range
        df_C['C_N-N-1'] = s1 
        df_C['C_N-N'] = s2 
        df_C['C_N-N+1'] = s3 
        df_C['Energy_mean'] = s4
        df_C['Energy2_mean'] = s5 
        df_C['Energy3_mean'] = s6     
        
        if save_C_tables is True:
            df.to_excel(filepath_og, index=False)
            df_C.to_excel(filepath_C, index=False)
        return [df, df_C]
    
    def calculate_P_values(self, df, df_C, save_C_tables=True, 
                           filepath_og="tmp_C_table_original_P-data.xlsx", 
                           filepath_C="tmp_C_table_interpolated_P-data.xlsx"):
        """
        Calculate P values, ratio and weights (\Pi values)from C data.
        
        Parameters
        ----------
        df: Data frame containing the C data in the columns C_N-N-1, C_N-N and C_N-N+1. 
        df_C: The dataframe containing the interpolated data at each N. 
        save_C_tables: Whether to save the output interpolated C data (default: False)
        filepath_og: File path to which the original C data with the new features calculated
            are stored (default: tmp_C_table_original.xlsx)  
        filepath_C: File path to which the interpolated C data with the new features 
            are stored (default: tmp_C_table_interpolated.xlsx)
        
        Returns
        -------
        df: The original data frame with the new features such as weights. 
        df_C: The dataframe containing the interpolated data with the new features at each N. 
        """
        
        #In this cell, we will manipulate the C values to get the weights at the reference pressure. 
        df_C['P_N-N-1'] = df_C['C_N-N-1'] / (df_C['C_N-N+1'] + df_C['C_N-N'] + df_C['C_N-N-1'])
        df_C['P_N-N+1'] = df_C['C_N-N+1'] / (df_C['C_N-N+1'] + df_C['C_N-N'] + df_C['C_N-N-1'])
        #First, we will take the log of the probabilities to ensure that they aren't too big or small. 
        df_C['log_P_N-N-1'] = scipy.log(df_C['P_N-N-1'])
        df_C['log_P_N-N+1'] = scipy.log(df_C['P_N-N+1'])
        #Next, we compute the ratios of \Pi values for the various loadings. This comes from detailed balance equations. 
        df_C['log_Pi_N+1/Pi_N'] = df_C['log_P_N-N+1']  -  df_C['log_P_N-N-1'].shift(-1) 
        df_C['log_Pi_N/Pi_N+1'] =   df_C['log_P_N-N-1'].shift(-1) -  df_C['log_P_N-N+1'] 
        #Or, amore elegant way. 
        df_C['log_Pi_N'] = 0+df_C['log_Pi_N+1/Pi_N'].cumsum().shift(1)
        df_C.loc[df_C['N']==0, 'log_Pi_N']=0
        
        #This won't work when we spacing!=1. Disregard this. 
        #In this cell, we will manipulate the C values to get the weights at the reference pressure for the reference C values. 
        df['P_N-N-1'] = df['C_N-N-1'] / (df['C_N-N+1'] + df['C_N-N'] + df['C_N-N-1'])
        df['P_N-N+1'] = df['C_N-N+1'] / (df['C_N-N+1'] + df['C_N-N'] + df['C_N-N-1'])
        #First, we will take the log of the probabilities to ensure that they aren't too big or small. 
        df['log_P_N-N-1'] = scipy.log(df['P_N-N-1'])
        df['log_P_N-N+1'] = scipy.log(df['P_N-N+1'])
        #Next, we compute the ratios of \Pi values for the various loadings. This comes from detailed balance equations. 
        df['log_Pi_N+1/Pi_N'] = df['log_P_N-N+1'] - df['log_P_N-N-1'].shift(-1) 
        df['log_Pi_N/Pi_N+1'] = df['log_P_N-N-1'].shift(-1) - df['log_P_N-N+1'] 
        #Or, amore elegant way. 
        df['log_Pi_N'] = 0+df['log_Pi_N+1/Pi_N'].cumsum().shift(1)
        df.loc[df['N']==0, 'log_Pi_N']=0

        if save_C_tables is True:
            df.to_excel(filepath_og, index=False)
            df_C.to_excel(filepath_C, index=False)
        return [df, df_C]
    
    def plot_C_values(self, df=None, df_C=None, fig=None, figure_title="", 
                    figure_description="C values", 
                    plot_data=True, plot_color="red", 
                    plot_label="", 
                    final_adjustment=True, save_image=False, 
                    image_filepath="tmp_iso_image.png", 
                    close_fig=False):
        """
        Plot C-related data using C values.

        Can be used in many ways using different combinations of parameters. This is 
        to initialize, plot or finish up by adjusting the axes and saving the figure. 
        All parameters are optional. 

        Parameters
        ----------
        df: The original dataframe to be plotted. 
        df_C: The interpolated dataframe to be plotted. 
        fig: Figure object, if present. 
        figure_title: Title of the figure. 
        figure_description: Description of the figures. 
        plot_data: Whether or not to plot the data.
        plot_color: Color of the data plotted. 
        plot_label: Label of the things plotted. 
        final_adjustment: After all the data has been plotted, the axes are readjusted 
            to include everuthing. 
        save_image: Whether or not to save figure to file. 
        image_filepath: File path to which figure is saved. 
        close_fig: Whether or not to close figure object. Only do this when plotting is finished. 
        
        Returns
        -------
        fig: Figure object.
        """
        
        if fig is None:
            #Figure for C values and spline fitting. 
            fig=plt.figure(figsize=(23, 6))
            fig.subplots_adjust(wspace=0.3)
            ax1=fig.add_subplot(131)
            ax2=fig.add_subplot(132)
            ax3=fig.add_subplot(133)
        
            #Also, initialize the axes with the necessary things. 
            ax1.set_xlabel('N')
            ax1.set_ylabel('C')
            ax2.set_xlabel('N')
            ax2.set_ylabel('C')
            ax3.set_xlabel('N')
            ax3.set_ylabel('C')

            ax1.set_title(figure_title)
            ax3.set_title(figure_description)
        else: #Figure exists. This means that it should also have the appropriate axes. 
            [ax1, ax2, ax3] = fig.get_axes()

        if plot_data is True:
            #N->N-1
            ax1.scatter(df['N'], df['C_N-N-1'], marker='o', color="white", edgecolor=plot_color, zorder=1.5) 
            ax1.plot(df_C["N"], df_C["C_N-N-1"], color=plot_color, linewidth=4, label=plot_label, zorder=1)
            
            #N->N
            ax2.scatter(df['N'], df['C_N-N'], marker='o', color="white", edgecolor=plot_color, zorder=1.5)
            ax2.plot(df_C["N"], df_C["C_N-N"], color=plot_color, linewidth=4, label=plot_label, zorder=1)
            
            #N->N+1
            ax3.scatter(df['N'], df['C_N-N+1'], marker='o', color="white", edgecolor=plot_color, zorder=1.5)
            ax3.plot(df_C["N"], df_C["C_N-N+1"], color=plot_color, linewidth=4, label=plot_label, zorder=1)
 
        if final_adjustment is True: #This means that this is the last plot in the series,and that we just want to 
                                    #align everything properly
            ax1.set_yscale('log')
            #ax1.set_ylim((1e-4, 1e5 ))
            #ax1.set_ylim((0 , 50 ))
            #ax1.legend(loc='lower left', bbox_to_anchor=(1.01, 0.01))
            #ax1.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.01))

            ax2.set_yscale('log')
            #ax.set_ylim((1e4, 1e8 ))
            #ax.legend(loc='lower left', bbox_to_anchor=(1.01, 0.01))
            ax2.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.01))

            ax3.set_yscale('log')
            #ax.set_ylim((1e-5, 1e7 ))
            #ax.legend(loc='lower left', bbox_to_anchor=(1.01, 0.01))
            #ax.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.01))

        if save_image is True:
            fig.savefig(image_filepath, dpi=300, bbox_inches='tight' )
        if close_fig is True:
            plt.close(fig)
            return None
        else:
            return fig

    def plot_P_values(self, df=None, df_C=None, fig=None, figure_title="", 
                    figure_description="P and Ratio", 
                    plot_data=True, plot_color="red", 
                    plot_label="", 
                    final_adjustment=True,save_image=True, 
                    image_filepath="tmp_iso_image.png",
                    close_fig=False):
        """
        Plot P data and ratio values.

        Can be used in many ways using different combinations of parameters. This is 
        to initialize, plot or finish up by adjusting the axes and saving the figure. 
        All parameters are optional. 

        Parameters
        ----------
        df: The original dataframe to be plotted. 
        df_C: The interpolated dataframe to be plotted. 
        fig: Figure object, if present. 
        figure_title: Title of the figure. 
        figure_description: Description of the figures. 
        plot_data: Whether or not to plot the data.
        plot_color: Color of the data plotted. 
        plot_label: Label of the things plotted. 
        final_adjustment: After all the data has been plotted, the axes are readjusted 
            to include everuthing. 
        save_image: Whether or not to save figure to file. 
        image_filepath: File path to which figure is saved. 
        close_fig: Whether or not to close figure object. Only do this when plotting is finished. 
        
        Returns
        -------
        fig: Figure object.
        """
        
        if fig is None:
            #Figure for C values and spline fitting. 
            fig=plt.figure(figsize=(23, 6))
            fig.subplots_adjust(wspace=0.3)
            ax1=fig.add_subplot(131)
            ax2=fig.add_subplot(132)
            ax3=fig.add_subplot(133)
        
            #Also, initialize the axes with the necessary things. 
            ax1.set_xlabel('N')
            ax1.set_ylabel('P_N_N-1')
            ax2.set_xlabel('N')
            ax2.set_ylabel('P_N_N+1')
            ax3.set_xlabel('N')
            ax3.set_ylabel('Ratio')

            ax1.set_title(figure_title)
            ax3.set_title(figure_description)
        else: #Figure exists. This means that it should also have the appropriate axes. 
            [ax1, ax2, ax3] = fig.get_axes()
     
        if plot_data is True:
            #N->N-1
            ax1.scatter(df['N'], df['C_N-N-1']/(df['C_N-N-1'] +df['C_N-N'] + df['C_N-N+1'] ), marker='o', color="white", edgecolor=plot_color, zorder=1.5)
            ax1.plot(df_C['N'], df_C['C_N-N-1']/(df_C['C_N-N-1'] + df_C['C_N-N'] + df_C['C_N-N+1']), color=plot_color, linewidth=4, label=plot_label, zorder=1)
            
            #N->N
            ax2.scatter(df['N'], df['C_N-N+1']/(df['C_N-N-1']+df['C_N-N']+df['C_N-N+1']), marker="o", color="white", edgecolor=plot_color, zorder=1.5)
            ax2.plot(df_C['N'], df_C['C_N-N+1']/(df_C['C_N-N-1']+df_C['C_N-N']+df_C['C_N-N+1']), color=plot_color, linewidth=4, label=plot_label, zorder=1)

            #N->N-1      
            ax3.scatter(df['N'], df['C_N-N-1']/df['C_N-N+1'], marker="o", color="white", edgecolor=plot_color, zorder=1.5)
            ax3.plot(df_C['N'], df_C['C_N-N-1']/df_C['C_N-N+1'], color=plot_color, linewidth=4, label=plot_label, zorder=1)

        if final_adjustment is True: #This means that this is the last plot in the series,and that we just want to 
                                    #align everything properly
            ax1.set_yscale('log')
            #ax1.set_ylim((1e-4, 1e5 ))
            #ax1.set_ylim((0 , 50 ))
            #ax1.legend(loc='lower left', bbox_to_anchor=(1.01, 0.01))
            #ax1.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.01))

            ax2.set_yscale('log')
            #ax2.set_ylim((1e4, 1e8 ))
            #ax2.legend(loc='lower left', bbox_to_anchor=(1.01, 0.01))
            ax2.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.01))

            ax3.set_yscale('log')
            #ax3.set_ylim((1e4, 1e8 ))
            #ax3.set_xlim((0, 1000))
            #ax3.legend(loc='lower left', bbox_to_anchor=(1.01, 0.01))
            #ax3.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.01))

        if save_image is True:
            fig.savefig(image_filepath, dpi=300, bbox_inches='tight' )
        if close_fig is True:
            plt.close(fig)
            return None
        else:
            return fig

    def plot_energy_values(self, df=None, df_C=None, 
                        fig=None, figure_title="", 
                        figure_description="Energy", 
                        plot_data=True, plot_color="red", 
                        plot_label="", 
                        final_adjustment=True,
                        save_image=False, 
                        image_filepath="tmp_iso_energy_image.png",
                        close_fig=False):
        """
        Plot energy and its moments using from the values.

        Can be used in many ways using different combinations of parameters. This is 
        to initialize, plot or finish up by adjusting the axes and saving the figure. 
        All parameters are optional. 

        Parameters
        ----------
        df: The original dataframe to be plotted. 
        df_C: The interpolated dataframe to be plotted. 
        fig: Figure object, if present. 
        figure_title: Title of the figure. 
        figure_description: Description of the figures. 
        plot_data: Whether or not to plot the data.
        plot_color: Color of the data plotted. 
        plot_label: Label of the things plotted. 
        final_adjustment: After all the data has been plotted, the axes are readjusted 
            to include everuthing. 
        save_image: Whether or not to save figure to file. 
        image_filepath: File path to which figure is saved. 
        close_fig: Whether or not to close figure object. Only do this when plotting is finished. 
        
        Returns
        -------
        fig: Figure object.
        """

        if fig is None:
            #Figure for C values and spline fitting. 
            fig=plt.figure(figsize=(23, 6))
            fig.subplots_adjust(wspace=0.3)
            ax1=fig.add_subplot(131)
            ax2=fig.add_subplot(132)
            ax3=fig.add_subplot(133)
        
            #Also, initialize the axes with the necessary things. 
            ax1.set_xlabel('N')
            ax1.set_ylabel('Energy mean')
            ax2.set_xlabel('N')
            ax2.set_ylabel('Energy2 mean')
            ax3.set_xlabel('N')
            ax3.set_ylabel('Energy3 mean')

            ax1.set_title(figure_title)
            ax3.set_title(figure_description)
        else: #Figure exists. This means that it should also have the appropriate axes. 
            [ax1, ax2, ax3] = fig.get_axes()

        if plot_data is True:
            ax1.scatter(df['N'], df['Energy_mean'], marker="o", color="white", edgecolor=plot_color, zorder=1.5)
            ax1.plot(df_C['N'], df_C['Energy_mean'], color=plot_color, linewidth=4, label=plot_label, zorder=1)
            
            ax2.scatter(df['N'], df['Energy2_mean'], marker="o", color="white", edgecolor=plot_color, zorder=1.5)
            ax2.plot(df_C['N'], df_C['Energy2_mean'], color=plot_color, linewidth=4, label=plot_label, zorder=1)
                    
            ax3.scatter(df['N'], df['Energy3_mean'], marker="o", color="white", edgecolor=plot_color, zorder=1.5)
            ax3.plot(df_C['N'], df_C['Energy3_mean'], color=plot_color, linewidth=4, label=plot_label, zorder=1)
        
        if final_adjustment is True: #This means that this is the last plot in the series,and that we just want to 
                                    #align everything properly
            #ax1.set_yscale('log')
            #ax1.set_ylim((1e-4, 1e5 ))
            #ax1.set_ylim((0 , 50 ))
            #ax1.legend(loc='lower left', bbox_to_anchor=(1.01, 0.01))
            #ax1.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.01))

            #ax2.set_yscale('log')
            #ax2.set_ylim((1e4, 1e8 ))
            #ax2.legend(loc='lower left', bbox_to_anchor=(1.01, 0.01))
            ax2.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.01))

            #ax3.set_yscale('log')
            #ax3.set_ylim((1e4, 1e8 ))
            #ax3.set_xlim((0, 1000))
            #ax3.legend(loc='lower left', bbox_to_anchor=(1.01, 0.01))
            #ax3.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.01))

        if save_image is True:
            fig.savefig(image_filepath, dpi=300, bbox_inches='tight' )
        if close_fig is True:
            plt.close(fig)
            return None
        else:
            return fig

    def reweight_iso(self, df, p, N_col='N', log_pi_N_muref_col='log_Pi_N', 
                    p_ref=4.11e9, T=298, drop_cols="Yes", T_ref = 298):
        """
        Reweight the log of the weights of the isotherms according to the formula. 
        
        Parameters:
        -----------
        df: The dataframe that should have \Pi values at the reference temperature and pressure.
        p: The pressures at which reweighting is to be performed. 
        T: Temperature list at which to reweight (default: 298)
        p_ref: Reference pressure (fugacity) at which simulation was done (default: 4.11e9 Pa)
        T_ref: Reference temperature at which simulation is carried out (default: 298).
        N_col: The column to be used for the N values (default: N).
        log_pi_N_muref_col: The column to be used for the log Pi at the reference 
            pressure (default'log_Pi_N'). 
        drop_cols: Whether or not to drop certain extra columns which may not be necessary 
            for calculations (default: "Yes")
        
        Returns
        -------
        df: Returns dataframe with the reweighted data at the required T and P conditions. 
        """
        #print ("Calculating for different T")
        #R = 8.314 #J/molK
        betadash = 1 / T # 1/K
        beta = 1 / T_ref # 1/K
        #dlnQdbeta = -df['Energy_mean']
        #d2lnQdbeta2 = df['Energy2_mean'] - df['Energy_mean']**2
        #d3lnQdbeta3 = - df['Energy3_mean'] - 2 * df['Energy_mean']**3 + 3 * df['Energy2_mean'] * df['Energy_mean']
        #lnQcdash_by_Qc = dlnQdbeta * (betadash - beta) +  1 / math.factorial(2) * d2lnQdbeta2 * (betadash - beta)**2 + 1 / math.factorial(3) * d3lnQdbeta3 * (betadash - beta)**3
        #lnQcdash_by_Qc = dlnQdbeta * (betadash - beta)  

        dlnQdbeta = -df['Energy_mean']
        #d2lnQdbeta2 = df['Energy2_mean'] - df['Energy_mean']**2
        #d3lnQdbeta3 = - df['Energy3_mean'] - 2 * df['Energy_mean']**3 + 3 * df['Energy2_mean'] * df['Energy_mean']
        lnQcdash_by_Qc = dlnQdbeta * (betadash - beta) #+  1 / math.factorial(2) * d2lnQdbeta2 * (betadash - beta)**2 + 1 / math.factorial(3) * d3lnQdbeta3 * (betadash - beta)**3

        df['log_pi_N_%.2e'%p] = df[log_pi_N_muref_col] + df[N_col] * scipy.log(p/p_ref) + lnQcdash_by_Qc
        df['scaled_log_pi_N_%.2e'%p]= df['log_pi_N_%.2e'%p] - df['log_pi_N_%.2e'%p].max() + 10 #The max log value we will have is 10.
        df['Pi_N_%.2e'%p] = scipy.exp(df['scaled_log_pi_N_%.2e'%p])
        df.loc['Total', 'Pi_N_%.2e'%p] = df['Pi_N_%.2e'%p].sum(axis=0)
        df['N_weights_%.2e'%p] = df[N_col] * df['Pi_N_%.2e'%p]
        df.loc['Total', 'N_weights_%.2e'%p] = df['N_weights_%.2e'%p].sum() / df.loc['Total', 'Pi_N_%.2e'%p]
        if drop_cols=="Yes": #This is to drop cols so that the dataframe becomes less cluttered.
            df.drop(['log_pi_N_%.2e'%p, 'scaled_log_pi_N_%.2e'%p], axis=1, inplace=True)                               
        return df
    
    def get_isotherm_data(self, df_C, temp_list=[298], pressure_list=[918], T_ref=298, P_ref=918, interp_method="regression", filename_root="rd2", 
                        make_weight_figure=False, save_weight_figure=False, weight_figure_filepath="tmp_weight_figure.png",
                        save_isotherm_file=True, isotherm_filepath='tmp_iso_reweighted.txt', 
                        save_full_C_data_file=True, C_data_filepath='tmp_C_data.xlsx'):
        """
        Compute isotherms from Pi data.
        
        **Make sure to set the "P_ref" and "interp method" same as that seen in the "interpolate_C_data" function manually. 
        
        Parameters
        --------------
        df_C: The dataframe that should have \Pi values at the reference temperature and pressure.
        temp_list: Temperature list at which to reweight (default: 298)
        pressure_list: The pressures at which reweighting is to be performed. 
        T_ref: Reference temperature at which simulation is carried out (default: 298).
        p_ref: Reference pressure (fugacity) at which simulation was done (default: 4.11e9 Pa)
        interp_method: Method via which interpolation is performed. Choose between "interpolation"
            or "regression" (default: regression). 
        filename_root: The root of the file name for identification.
        make_weight_figure: Whether or not to make a figure of weights at each required 
            temperature or pressure. 
        save_weight_figure: Save the weight figure if created. 
        weight_figure_filepath: File path to which to save figure is figure is to be saved. 
        save_isotherm_file: Whether to save isotherm file. 
        isotherm_filepath: File path of the isotherm file. 
        save_full_C_data_file: Whether the C data file with weights at all desired T and P is 
            to be saved. 
        C_data_filepath: Filepath of the C data if it is to be saved. 
        
        Returns
        --------------
        df_iso: Isotherm dataframe.
        df_C: Dataframe with weights at all desired T and P conditions. 
        """
        for T in temp_list:
            #Here, we will extract the loading values for the pressures. 
            df_iso= pd.DataFrame(columns=['Pressure', 'NLoading'])
            
            for index, p in enumerate ( pressure_list ):
                #print ("Calculating for pressure = %.2e, %d"%(p, T))
                df_iso.loc[index, 'Pressure'] = p
                df_C = self.reweight_iso(df_C, p, drop_cols="No", p_ref = P_ref, T=T, T_ref=T_ref)
                df_iso.loc[index, 'NLoading'] = df_C.loc['Total', 'N_weights_%.2e'%p]
                                
                if make_weight_figure is True: #For diagnostic purposes.
                    fig=plt.figure()
                    ax=fig.add_subplot(111)
                    ax.scatter(df_C['N'], df_C['scaled_log_pi_N_%.2e'%p] - df_C['scaled_log_pi_N_%.2e'%p].max() )
                    ax.set_title('Pressure= %.2e'%p)
                    ax.set_xlabel('N')
                    ax.set_ylabel('scaled_log_pi_N')
                    #ax.set_ylim((-350, 10))
                    ax.set_ylim((-50, 10))
                    if save_weight_figure is True:
                        fig.savefig(weight_figure_filepath, dpi=300, bbox_inches="tight")
                                
            if save_isotherm_file is True:
                df_iso.to_csv(isotherm_filepath, index=False)
            if save_full_C_data_file is True:
                df_C.to_excel(C_data_filepath, index=False)

        return [df_iso, df_C]
    
    def plot_isotherm_figure(self, df_iso=None, fig=None, ax=None,
                        figure_description="Isotherms", 
                        axes_labels=["Fugacity (Pa)", "Loading (molecules/framework)"],
                        plot_data=True, 
                        x_value="Pressure", 
                        y_value="NLoading",
                        plot_color="red", 
                        plot_label="", final_adjustment=True, 
                        xlim_values=None,
                        save_isotherm_figure=True, 
                        isotherm_figure_filepath="tmp_iso_figure.png",
                        close_fig=False):
        """
        Plot isotherm figure.
        
        Can be used in many ways using different combinations of parameters. This is 
        to initialize, plot or finish up by adjusting the axes and saving the figure. 
        All parameters are optional. 

        Even axes can be directly worked upon if "fig" parameter is set to None and 
        "ax" is parsed.   

        Parameters
        ----------
        df_iso: Isotherm dataframe.
        fig: Figure object in which data can be plotted. 
        ax: Axes objects which can be acted upon. 
        figure_description: Description of the figures. 
        axes_labels: Labels for the axes.
        plot_data: Whether or not to plot the data.
        x_value: The value to be plotted on the x axis.
        y_value: The value to be plotted on the y axis.
        plot_color: Color of the data plotted. 
        plot_label: Label of the things plotted. 
        final_adjustment: After all the data has been plotted, the axes are readjusted 
            to include everuthing. 
        xlim_values: X lim values (Optional, default: None).
        save_isotherm_figure: Whether or not to save figure to file. 
        isotherm_figure_filepath: File path to which figure is saved. 
        close_fig: Whether or not to close figure object. Only do this when plotting is finished. 
        
        Returns
        -------
        fig or ax: Figure or axis object with actions as specified. 
        """
        #if fig is None:    
        if ax is None and fig is None:  #Neither axis nor object provided. make figure.
            fig=plt.figure()
            ax=fig.add_subplot(111)
            ax.set_xlabel(axes_labels[0])
            ax.set_ylabel(axes_labels[1])
            ax.set_title(figure_description)
        elif ax is None and fig is not None: #Figure already provided. extract axis. 
            [ax]=fig.get_axes()
        elif ax is not None and fig is None:#the other possibility is that ax is 
                            #present and figure is None. In that case, we just move on. 
            pass

        if plot_data is True:
            ax.plot(df_iso[x_value], df_iso[y_value], color=plot_color, label=plot_label, linewidth=4)

        if final_adjustment is True:
            ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))

        if xlim_values is not None:
            ax.set_xlim(xlim_values)

        if fig is not None: #Figure object hasn't been instantiated. 
            if save_isotherm_figure is True:
                fig.savefig(isotherm_figure_filepath, dpi=300, bbox_inches="tight")
            if close_fig is True:
                plt.close(fig)
            return fig
        else: #figure is None, that means we only want to modify the axis.
            return ax
