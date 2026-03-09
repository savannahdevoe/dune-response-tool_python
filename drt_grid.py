# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 12:41:32 2025

@author: sdevo
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d

def drt_grid(scenario):
    #scenario = drt_grid(scenario)
    # drt_erosion: code to set up grids for input to the dune response tool
    
    # # Build Profile 
    # # assumed geometric and gridding parameters for simplicity
    # scenario['grids'] = {} # already done in drt_search_morphology
    # scenario['grids']['morphometrics'] = {} # already done in drt_search_morphology
    scenario['grids']['morphometrics']['dunewidth'] = 10 
    scenario['grids']['morphometrics']['dx'] = 0.1
    scenario['grids']['morphometrics']['smoothlev'] = 10 # divide by 5 to get gaussian filter sigma to equal the window size in original matlab DRT
    scenario['grids']['morphometrics']['shore'] = 0

    # set up profile segments
    numdune = np.int16(np.round(scenario['grids']['morphometrics']['dunewidth']/scenario['grids']['morphometrics']['dx']))
    Zbackdune = (scenario['grids']['morphometrics']['dhigh']-2)*np.ones(np.int16(np.ceil(20/scenario['grids']['morphometrics']['dx'])))
    Zbackduneslope = np.arange(scenario['grids']['morphometrics']['dhigh']-0.5,scenario['grids']['morphometrics']['dhigh'],scenario['grids']['morphometrics']['duneslope']*scenario['grids']['morphometrics']['dx'])
    Zdune1 = scenario['grids']['morphometrics']['dhigh']*np.ones((numdune))
    Zdune2 = np.arange(scenario['grids']['morphometrics']['dhigh'],scenario['grids']['morphometrics']['dtoe'],-(scenario['grids']['morphometrics']['duneslope'])*scenario['grids']['morphometrics']['dx'])
    Zbeach = np.arange(scenario['grids']['morphometrics']['dtoe'],scenario['grids']['morphometrics']['shore'],-(scenario['grids']['morphometrics']['backshoreslope'])*scenario['grids']['morphometrics']['dx'])                  
    
    # combine all segments
    Zall =np.concatenate((Zbackduneslope[:], Zdune1[:], Zdune2[:], Zbeach[:]))
    Zall = np.flipud(Zall)

    #Smooth Data
    Zall = gaussian_filter1d(Zall, sigma=scenario['grids']['morphometrics']['smoothlev'])

    # Grid Setup
    XGrid = -(np.arange(0,len(Zall),1))*np.abs(scenario['grids']['morphometrics']['dx'])
    XGrid = np.flipud(XGrid).T
    
    # Output Data
    scenario['grids']['XGrid'] = XGrid.T
    scenario['grids']['ZGrid'] = Zall.T
    
    return scenario