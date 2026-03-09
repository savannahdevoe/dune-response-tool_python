# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 09:10:24 2025

@author: sdevo
"""
import numpy as np
# import math

def drt_accretion(scenario):
    # drt_accretion: code to run a simple aeolian sediment transport model
    # for calculating wind blown sediment fluxes into coastal dunes
    #
    # Required Inputs: 'scenario' dict variable with the 'grids', 
    # 'erosion', and 'env' dict variables
    #
    # Outputs: scenario.accretion
    
    # set up aeolian model
    twl = scenario['erosion']['TWL']
    xprof = scenario['grids']['XGrid']
    zprof = scenario['grids']['ZGrid']
    dtoe = scenario['grids']['morphometrics']['dtoe']
    dhigh = scenario['grids']['morphometrics']['dhigh']
    imax = np.argmax(zprof)
    xtoe = np.interp(dtoe,zprof[0:imax+1],xprof[0:imax+1])
    u_w = scenario['env']['winds']['windSpeed'] # [m/s]
    windDir = scenario['env']['winds']['windDirection']-scenario['grids']['morphometrics']['azimuth']
    
    # Utilizing Kawamura (1951) approach for wind-driven sediment fluxes
    D50 = scenario['models']['d50']/1000 # grain size
    K = 0.4 # von Karman constant
    z = 10 # assumed elevation of wind measurements
    zo = 2*D50/30 # roughness length scale (?)
    pa = 1.225 # air density [kg/m^3]
    ps = 2650 # sediment density [kg/m^3]
    g = 9.81 # gravity [m/s^2]
    u_star = u_w*K/np.log(z/zo) # shear velocity [m/s]
    C = 1.87 # for typical sands
    M = 0 # assumed moisture content
    Ck = scenario['models']['AeolianTransportCoefficient'] # model coefficient now user input
    ustar_thresh = 0.1* np.sqrt(g*D50*(ps/pa)*(1+C*M)) # threshold shear velocity [m/s]
    Q = Ck*(pa/g)*(u_star-ustar_thresh)*(u_star+ustar_thresh)*(u_star+ustar_thresh) # potential volumetric transport rate [m^3/m/s]
    
    # modify transport by the fetch effect per Delgado-Fernandez, Geomorphology (2011)
    # critical fetch length:
    Fc = 4.38*u_w - 8.23
    
    # determine the beach width based on the total water level
    maxval, imax = np.nanmax(zprof), np.argmax(zprof)
    minval = np.nanmin(zprof)
    beachwidth = np.full(np.size(twl),np.nan)
    for itime in range(0,np.size(twl)):
        if twl[itime]< maxval and twl[itime] > minval: # if twl is within the profile
            xwl = np.interp(twl[itime],zprof[0:imax+1],xprof[0:imax+1])
            beachwidth[itime] = xtoe-xwl
        elif twl[itime] <= minval: # if twl is lower than the profile goes, set beach width to be wide:
            beachwidth[itime] = 100 # pick some high number if exceeded
        else:
            beachwidth[itime] = 0 
    
    beachwidth[np.where(beachwidth<0)] = 0 
    
    # determine the fetch for the specific conditions:
    F = beachwidth*[np.tan(ang) for ang in [np.radians(ang) for ang in np.abs(windDir)]]
    F[np.where(F<0)] = 0
    F[np.isinf(F)] = 1000
    Qtot = Q # intialize variable
    for idx in range(0,np.size(Qtot)):
        if F[idx] < Fc[idx]:
            Qtot[idx] = Q[idx]*np.sin(((np.pi/2)*F[idx]/Fc[idx]))
    Qtot[np.where((Qtot<0))] = 0
    
    # lastly, deal with flux to dunes based on angle to get flux to dune
    Qtot = Qtot*[np.cos(ang) for ang in [np.radians(ang) for ang in np.abs(windDir)]]
    Qtot[np.where(Qtot<0)] = 0
    
    # convert to a volume flux (initially in kg/m/s)
    por = 0.4 # assumed porosity
    Qtot_m3m_dt = (Qtot/ps)*(scenario['timing']['dt']*60*60)/(1-por)
    
    # convert outputs
    scenario['accretion'] = {}
    scenario['accretion']['dV'] = Qtot_m3m_dt
    
    return scenario