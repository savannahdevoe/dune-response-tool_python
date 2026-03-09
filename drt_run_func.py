# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 09:25:15 2025

@author: sdevo
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 13:11:47 2025

@author: sdevo
"""
# import required functions and packages
import drt_env
import drt_grid
import drt_erosion
import drt_accretion
import drt_plotting

import pandas as pd
import datetime
from zoneinfo import ZoneInfo
import pickle
    
def drt_run_func(scenario,plotting=True,savefile=False):
    # # # End Date (Automatically Calculated) # # #
    scenario['timing']['end_date'] = scenario['timing']['start_date'] + datetime.timedelta(days=scenario['timing']['Durationdays'])
    
    # # # Set up grids based on dune and beach morphology conditions: # # #
    # look for dune morpho from USGS data based on closest lat/lon:
    scenario = drt_env.drt_search_morphology(scenario, morph_file=f"{scenario['code_direc']}\dependencies\drt_morphology.xlsx")
    # make model x and z grids:
    scenario = drt_grid.drt_grid(scenario)
    
    # # # Pull Environmental Variables (Wind, Waves, Tides) # # #
    scenario = drt_env.drt_env(scenario)
    
    # # # Run Erosion Module # # #
    scenario = drt_erosion.drt_erosion(scenario)
    
    # # # Run Accretion Module # # #
    scenario = drt_accretion.drt_accretion(scenario)
    
    if plotting:
        # # # Plot the Outputs # # #
        drt_plotting.drt_plotting(scenario)
        
    if savefile:
        # # # Save an Output File # # #
        with open(f"DRT_{scenario['timing']['start_date'].strftime('%Y%m%d')}-{scenario['timing']['end_date'].strftime('%Y%m%d')}_{scenario['location']['lat']}_{scenario['location']['lon']}.pkl",'wb') as fileout:
            pickle.dump(scenario,fileout)
    
    return scenario