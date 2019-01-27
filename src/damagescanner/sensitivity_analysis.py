"""DamageScanner - a directe damage assessment toolkit

Sensitivity analysis using the damagescanner
"""

import os
import pandas
import numpy

from SALib.sample import morris
import SALib.analyze.morris

from tqdm import tqdm
import matplotlib.pyplot as plt

import rasterio

from damagescanner.core import RasterScanner

def run(landuse_map,inun_map,curve_path,maxdam_path,save=False,**kwargs):
    """
    Perform sensitivity analysis of the results. Currently only works with the RasterScanner.
    
     Arguments:
        *landuse_map* : Shapefile, Pandas DataFrame or Geopandas GeoDataFrame 
        with land-use information of the area.
     
        *inun_map* : GeoTiff with inundation depth per grid cell. Make sure 
        that the unit of the inundation map corresponds with the unit of the 
        first column of the curves file. 
     
        *curve_path* : File with the stage-damage curves of the different 
        land-use classes. Can also be a pandas DataFrame (but not a numpy Array).
     
        *maxdam_path* : File with the maximum damages per land-use class 
        (in euro/m2). Can also be a pandas DataFrame (but not a numpy Array).
        
    Returns:
        *Figure* : Functions returns a figure with the spread in losses and 
        a spyder plot with the influence of each parameter.

    """
    
    # Set up problem for sensitivity analysis
    problem = {
          'num_vars': 3,
          'names': ['inundation', 'maxdam','curves'],
          'bounds': [[0.5,1.5],[0.5,1.5],[1,5]]}
    
    # And create parameter values
    param_values = morris.sample(problem, 10, num_levels=5, grid_jump=2,local_optimization =True)
    
    # Run analysis with the specified parameter values
    collect_losses = []
    for param_set in tqdm(param_values,total=len(param_values)):
        with rasterio.open(inun_map) as src:
            inundation = src.read()[0,:,:]*param_set[0]
        inundation[inundation>10] = 0
        inundation[inundation<0] = 0
        
        maxdam = pandas.read_csv(maxdam_path,skiprows=1).values*param_set[1]
        
        curves = pandas.read_csv(curve_path).values
        curves[:,1:] = curves[:,1:]/param_set[2]
        
        loss_df = RasterScanner(landuse_map,inundation,curves,maxdam,save=False,
                                                         in_millions=True)
        collect_losses.append(loss_df[0].sum().values[0])
        
    # Obtain results
    Si = SALib.analyze.morris.analyze(problem, numpy.array(param_values), 
                                          numpy.array(collect_losses),
                                     print_to_console=False, grid_jump=2, num_levels=5)
    
    # And plot the results
    fig = plt.figure(figsize=(10, 5))
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122, projection='polar')
    
    color_scheme_map = ['#264653','#2A9D8F']
        
    # Plot boxplot with spread of total losses
    ax1.boxplot(numpy.array(collect_losses),showfliers=False, patch_artist=True,
                boxprops=dict(facecolor=color_scheme_map[1], edgecolor='black'),
                capprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                flierprops=dict(color=color_scheme_map[1], markeredgecolor='black'),
                medianprops=dict(color=color_scheme_map[1]),
                )
    
    ax1.set_title('Spread of total losses',fontsize=18,fontweight='black', y=1.08)
    
    # Plot spyder plot with influence of different variables
    risk_sens = pandas.DataFrame.from_dict(Si)
    risk_sens['rel'] = abs(risk_sens['mu'])/abs(risk_sens['mu']).sum()*100
    risk_sens['mu'] = risk_sens['mu']
    risk_sens = risk_sens.groupby('names').sum()
    risk_sens = risk_sens.T
    risk_sens.columns = ['inundation', 'maxdam', 'curve']      
    
    stats=risk_sens.loc['rel',numpy.array(risk_sens.columns)].values
    
    angles=numpy.linspace(0, 2*numpy.pi, len(numpy.array(risk_sens.columns)), endpoint=False)
    
    # close the plot
    stats=numpy.concatenate((stats,[stats[0]]))
    angles=numpy.concatenate((angles,[angles[0]]))
    
    ax2.plot(angles, stats, 'o-', linewidth=2,color=color_scheme_map[0])
    ax2.set_ylim([0, 100])   
    ax2.fill(angles, stats, alpha=0.25,color=color_scheme_map[0])
    ax2.set_thetagrids(angles * 180/numpy.pi, numpy.array(risk_sens.columns))
    ax2.tick_params(axis='x',labelsize=14,labelcolor='black',color='black',) # pad=12
    
    ax2.set_title('Influence of variables',fontsize=18,fontweight='black', y=1.08)
    
    fig.tight_layout()
    
    # Save results if desired
    if save:
        output_path = kwargs['output_path']
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        fig.savefig(os.path.join(output_path,'sensitivity_analysis.png'),dpi=300)
    