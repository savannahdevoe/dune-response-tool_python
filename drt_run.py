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

# Main Run Segment 

# # # PREALLOCATE # # #
scenario = {}
scenario['location'] = {}
scenario['timing'] = {}
scenario['models'] = {}
scenario['env'] = {}

# # # USER INPUTS: # # #
scenario['code_direc'] = 'F:\PyPath\duneResponseTool' # python directory for DRT functions and dependables
scenario['location']['lat'] = 38.29#36.19 # this is FRF duneLidar loc. #39.44 # this is DRT example in NJ # scenario latitude (deg. N)
scenario['location']['lon'] = -75.11#-75.75 # this is FRF duneLidar loc. #-74.35 # this is DRT example in NJ # scenario longitude (deg. E.. note, negative = W)
scenario['type'] = 'hindcast' # 'forecast' or 'hindcast'?
scenario['timing']['start_date'] = datetime.datetime(1998,1,15)#datetime.datetime(2015,9,20)#,tzinfo=ZoneInfo('US/Eastern')) # enter start date
Durationdays = 30 # how many days do you want the simulation to run?
scenario['timing']['dt'] = 1 # model dt (hours)
scenario['models']['d50'] = 0.30#0.35 # median grain size diameter (mm)

# # # Model Coefficients (can be changed but these are defaults) # # #
scenario['models']['WaveRunupFactor'] = 1.26
scenario['models']['DuneSlopeTrajectory'] = 0.54
scenario['models']['DuneErodibility'] = 0.0025
scenario['models']['AeolianTransportCoefficient'] = 2.78
             
# # # End Date (Automatically Calculated) # # #
scenario['timing']['end_date'] = scenario['timing']['start_date'] + datetime.timedelta(days=Durationdays)

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

# # # Plot the Outputs # # #
drt_plotting.drt_plotting(scenario)

# # # Save an Output File # # #
with open(f"DRT_{scenario['timing']['start_date'].strftime('%Y%m%d')}-{scenario['timing']['end_date'].strftime('%Y%m%d')}_{scenario['location']['lat']}_{scenario['location']['lon']}.pkl",'wb') as fileout:
    pickle.dump(scenario,fileout)
    