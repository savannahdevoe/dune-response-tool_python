# -*- coding: utf-8 -*-
"""
Created on Mon Nov  3 16:05:54 2025

@author: sdevo
"""

import numpy as np
from scipy.stats import norm
import pandas as pd

def runPH12(xM,z,time,WL,Ho,Lo,T,Bo,zb,output_times, D50, WaveRunupFactor, DuneSlopeTrajectory, DuneErodibility):
    # PH12 = runPH12(xM,z,time,WL,Ho,Lo,T,Bo,zb,output_times, D50, WaveRunupFactor, DuneSlopeTrajectory, DuneErodibility)
    # Function to run the Palmsten and Holman (2012) model for dune retreat
    #   inputs:
    #       xM - x vector (with origin offshore)
    #       z - elevation profile (same size as xM)
    #       time - time vector wave data in seconds
    #       WL - water level time-series (tide + surge)
    #       Ho - deep water wave height time-series 
    #       Lo - deep water wave length time-series
    #       T - peak wave period time-series
    #       Bo - initial beach slope
    #       zb - initial dune toe elevation 

    # Set Model Coefficients
    d50 = D50/1000 # mm to m
    nsigma = 2 # in definition of R2, R16.. For R2, nsigma=2, R16 = nsigma=1; 
    g = 9.81 # gravity
    Kd = WaveRunupFactor # coefficient to account for higher runup on dune
    Cs = DuneErodibility 
    Ac = 1.34*10**(-3)
    bc = 3.19*10**(-4)
    Btfac = DuneSlopeTrajectory

    # Model Initialization
    val, st1 = np.min((z-zb)**2), np.argmin((z-zb)**2) # find grid point where initial dune toe is
    val2, st2 = np.min((z-zbd)**2), np.argmin((z-zbd)**2) # find grid point where initial back side of dune toe is
    Bt = Bo*Btfac # slope at which beta receeds. LEH04 = 1, PH11 = 0.54....
    
    zbT = np.concatenate((np.full((st1),np.nan), Bt*(xM[st1:]-xM[st1])+ zb))  # trajectory that dune toe recedes
    zbT = np.concatenate((np.full((st1),np.nan), Bt*(xM[st1:]-xM[st1])+ zb))  # trajectory that dune toe recedes

    ireplace = np.where(zbT>z)[0]
    zbT[ireplace] = z[ireplace]
    dt = np.diff(time[0:2]) # dt in seconds
  
    zb2 = np.full(len(WL)+1,np.nan)
    zb2[0] = zb
    zb = zb2
    
    # Main Program Loop
    output_num = 0

    xToe = np.full(len(WL)+1,np.nan)
    xBack = np.full(len(WL)+1,np.nan)
    V = np.full(len(WL),np.nan)
    Beta = np.full(len(WL),np.nan)
    etabar = np.full(len(WL),np.nan)
    sigma_s = np.full(len(WL),np.nan)
    zR = np.full(len(WL),np.nan)
    zRLEH = np.full(len(WL),np.nan)
    zTotal = np.full(len(WL),np.nan)
    p = np.full(len(WL),np.nan)
    Nc = np.full(len(WL),np.nan)
    dV = np.full(len(WL),np.nan)
    dVT = np.full(len(WL),np.nan)
    dx = np.full(len(WL),np.nan)
    dVResidual = np.full(len(WL),np.nan)
    Vover = np.full(len(WL),np.nan)
    zNew = []#np.full(len(WL),np.nan)
    actual_output_times = []#np.full(len(WL),np.nan)
  
    for tt in range(0,len(WL)):
        current_output_time = output_times[output_num]

        if tt==0:
            st = st1
            stt = st2
        else:
            st = st+ii
            stt = stt+ii
            
        # dune toe position
        xToe[tt] = xM[st]
        # back of dune toe position
        xBack = xm[stt]

        # # dune volume
        # V[tt] = np.sum(np.abs(np.diff(xM[0:2]))*(z[st:])) # measured in ref to z=0
        # Vc = np.cumsum(np.abs(np.diff(xM[0:2]))*(z[st:]-zbT[st:])) # cumulative volume above the dune trajectory
        V[tt] = np.sum(np.abs(np.diff(xM[st-1:]))*(z[st:])) # measured in ref to z=0
        Vc = np.cumsum(np.abs(np.diff(xM[st-1:]))*(z[st:]-zbT[st:])) # cumulative volume above the dune trajectory

        Vc = Vc - Vc[0]

        Beta[tt] = Bo # initial dummy guess.

        # stockdon for TWL
        etabar[tt] = 0.35*Beta[tt]*np.sqrt(Ho[tt]*Lo[tt]) # mean swash (setup)
        sigma_s[tt] = (np.sqrt(Ho[tt]*Lo[tt]*(0.563*(Beta[tt]**2)+0.0004))/2)*nsigma/2 # swash (IG + incident)
        zR[tt] = 1.1*(etabar[tt] + sigma_s[tt])
        # sigma_s2[tt] = np.sqrt(Ho[tt]*Lo[tt]*(0.563*(Beta[tt]**2)+0.0004))/2
        # zR[tt] = 1.1*(etabar[tt] + sigma_s[tt])

        zRLEH[tt] = 0.158*np.sqrt(Ho[tt]/1.416*Lo[tt]) # Larson et al. (2004) approximation
        zTotal[tt] = zR[tt]*Kd + WL[tt] # TWL is tidal + runup*runupfactor
        if zTotal[tt] >= np.max(z): # if TWL > max bed elevation
            zTotal[tt] = np.max(z)
        
        # R2 is 2% exceedance value for runup, calculated from CPDF of runup elevations
        # use mean value as WL, std is swash range
        p[tt] = 1 - norm.cdf(zb[tt],loc=etabar[tt]+WL[tt],scale=sigma_s[tt])
        Nc[tt] = p[tt]*(dt/T[tt])[0] # probability * dt/wave period is period of exposure to runup (proxy for # of collisions)
        if tt>0:
            # change in dune volume following Larson et al. (2004):
            dV[tt] = 4*Cs*(np.max(zTotal[tt]-zb[tt]))**2*Nc[tt]
            dVT[tt] = dV[tt] - dVResidual[tt-1] # volume change - residual volume left in the dune
        else:
            dVT[tt] = 4*Cs*(np.max(zTotal[tt]-zb[tt]))**2*Nc[tt]
            dV[tt] = 0
            
        if dVT[tt]<0: # if less than 0 vol. left in dune:
            ii=0
        else: # otherwise,
            val, ii = np.min((Vc-dVT[tt])**2), np.argmin((Vc-dVT[tt])**2) # find grid point where dune toe is
        dx[tt] = xM[st+ii]-xToe[tt] # change in x of dune toe position
        dVResidual[tt] = Vc[ii]-dVT[tt] # residual volume = dune volume - total volume change
        zb[tt+1] = Bt*dx[tt] + zb[tt] # trajectory that dune toe receeds.

        # overwash volume?
        xc = (xToe[tt]*np.tan(Bo) + xBack[tt]*np.tan(Bb))/(np.tan(Bo)+np.tan(Bb))
        Vover[tt] = dVResidual[tt] - 0.5*(xBack[tt]-xToe[tt])*(xc-xToe[tt])*np.tan(Bo)
        
        
        if time[tt] >= current_output_time:
            
            # zNew[output_num] = np.concatenate((z[0:st1], zbT[st1+1:st], z[st+1:])) # assumes vertical cliff face; probably needs a slope adjustmnet
            # actual_output_times[output_num] = time[tt]
            zNew.append(np.concatenate((z[0:st1+1], zbT[st1+1:st+1], z[st+1:]))) # assumes vertical cliff face; probably needs a slope adjustmnet
            actual_output_times.append(time[tt])
            output_num += 1

        # clean-up variables
        del Vc

    xToe[tt+1] = xToe[tt]+dx[tt]
    st = st+ii-1
    zbOut = zb[0:-1]
    xToe = xToe[0:-1]
    
    # dV[0] = 0
    # dVT[0] = 0
    # save a dict of model results
    PH12 = {}
    PH12['TWL'] = zTotal # final profile assuming a vertical cliff face (can under-estimate volumes of erosion)
    PH12['zNew'] = zNew # profiles at each time-step assuming a vertical cliff face (can under-estimate volumes of erosion)
    PH12['ztoe'] = zbOut # elevation of dune toe through time at each time step
    PH12['xtoe'] = xToe # cross-shore position of dune toe at each time step
    PH12['Nc'] = Nc # number of dune collisions at each time step
    PH12['dV'] = dV # volume of eroded sand at each time step
    PH12['dVResidual'] = dVResidual
    PH12['dVT'] = dVT
    PH12['zRunup'] = zR;
    PH12['times'] = time
    PH12['zmat_time'] = actual_output_times

    return PH12

def drt_erosion(scenario):
    # scenario = drt_erosion(scenario)
    # drt_erosion: code to run the Palmsten and Holman (2012) dune retreat
    # model
      
    # initialize model in PH12 conventions
    xM = scenario['grids']['XGrid']+np.max(np.abs(scenario['grids']['XGrid']))
    z = scenario['grids']['ZGrid']
    time = (scenario['timing']['times'] - pd.Timestamp("1970-01-01"))/pd.Timedelta("1s")
    time = time-time[0] # converted to seconds starting at t0 = 0
    WL = scenario['env']['tides']['wl']
    Ho = scenario['env']['waves']['Hs_25m']
    Lo = scenario['env']['waves']['L_25m']
    T = scenario['env']['waves']['Tp']
    Bo = scenario['grids']['morphometrics']['backshoreslope']
    zb = scenario['grids']['morphometrics']['dtoe']
    D50 = scenario['models']['d50']
    WaveRunupFactor = scenario['models']['WaveRunupFactor']
    DuneSlopeTrajectory = scenario['models']['DuneSlopeTrajectory']
    DuneErodibility = scenario['models']['DuneErodibility'] 

    # run model
    output_times = np.linspace(time[0],time[-1],10)    
    
    # output_times 
    PH12 = runPH12(xM,z,time,WL,Ho,Lo,T,Bo,zb, output_times, D50, WaveRunupFactor, DuneSlopeTrajectory, DuneErodibility)

    scenario['erosion'] = {}
    # store relevant model output and send back to main program
    scenario['erosion']['dV'] = PH12['dV']
    scenario['erosion']['Z'] = PH12['zNew']
    scenario['erosion']['TWL'] = PH12['TWL']
    scenario['erosion']['xToe'] = PH12['xtoe']
    scenario['erosion']['ztoe'] = PH12['ztoe']
    scenario['erosion']['Nc'] = PH12['Nc']
    scenario['erosion']['dVR'] = PH12['dVResidual']
    scenario['erosion']['dVT']= PH12['dVT']
    scenario['erosion']['times'] = PH12['times']
    scenario['erosion']['zmat_times'] = []
    for ii in range(len(PH12['zmat_time'])):
        temp = pd.Timedelta(value=PH12['zmat_time'][ii],unit='s')
        scenario['erosion']['zmat_times'].append(temp+scenario['timing']['times'][0])
    scenario['erosion']['zRunup'] = PH12['zRunup']

    return scenario