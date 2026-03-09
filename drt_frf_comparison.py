# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 10:45:47 2025

@author: sdevo
"""

import numpy as np
import datetime
from geoprocess import FRFcoord
import drt_run_func
import pickle
from download_FRFThredd_ncml import download_FRF_ncml
from concat_ncml_daterange import concat_ncml_daterange

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import gridspec

xrange = 42
yrange = np.arange(0,1000,25)

# # # PREALLOCATE # # #
scenario = {}
scenario['location'] = {}
scenario['timing'] = {}
scenario['models'] = {}
scenario['env'] = {}

# # # USER INPUTS: # # #
scenario['code_direc'] = 'F:\PyPath\duneResponseTool' # python directory for DRT functions and dependables
scenario['type'] = 'hindcast' # 'forecast' or 'hindcast'?
scenario['timing']['start_date'] = datetime.datetime(2015,9,20)#,tzinfo=ZoneInfo('US/Eastern')) # enter start date
scenario['timing']['Durationdays'] = 21 # how many days do you want the simulation to run?
scenario['timing']['dt'] = 1 # model dt (hours)
scenario['models']['d50'] = 0.35 # median grain size diameter (mm)

# # # Model Coefficients (can be changed but these are defaults) # # #
scenario['models']['WaveRunupFactor'] = 1.26
scenario['models']['DuneSlopeTrajectory'] = 0.54
scenario['models']['DuneErodibility'] = 0.0025
scenario['models']['AeolianTransportCoefficient'] = 2.78

cumDV_net = []
cumDXTOE = []
xTOE_start = []
xTOE_end = []
xgr = []
zgr = []
for ii, yy in enumerate(yrange):
    coords = FRFcoord(xrange,yy) # compute latitude and longitude from FRF coordinates
    
    scenario['location']['lat'] = coords['Lat'] # scenario latitude (deg. N)
    scenario['location']['lon'] = coords['Lon'] # scenario longitude (deg. E.. note, negative = W)

    scenario = drt_run_func(scenario,plotting=False)
    
    # if ii==0:
    #     fullZ = np.full((len(yrange),len(scenario['erosion']['Z'][0]),len(scenario['erosion']['Z'])),np.nan)
        
    cumDV_erosion = -np.cumsum(scenario['erosion']['dV'])[-1]
    cumDV_accretion = np.cumsum(scenario['accretion']['dV'])[-1]
    cumDV_net.append(cumDV_accretion + cumDV_erosion)
    xTOE_start.append(scenario['erosion']['xToe'][0]-np.max(np.abs(scenario['grids']['XGrid'])))
    xTOE_end.append(scenario['erosion']['xToe'][-1]-np.max(np.abs(scenario['grids']['XGrid'])))
    cumDXTOE.append(xTOE_end[-1]-xTOE_start[-1])
    # fullZ[ii,:,:] = np.asarray(scenario['erosion']['Z']).T
    xgr.append((scenario['grids']['XGrid']))
    zgr.append((scenario['erosion']['Z']))

cumDV_net = np.asarray(cumDV_net)
cumDXTOE = np.asarray(cumDXTOE)
xTOE_start = np.asarray(xTOE_start)
xTOE_end = np.asarray(xTOE_end)

fig = plt.figure(figsize=(16, 9))
gs = gridspec.GridSpec(5, 3, figure=fig, wspace=0.25, hspace=0.35)
r=[[2,0],[2,1],[2,2],[3,0],[3,1],[3,2],[4,0],[4,1],[4,2]]

for xx, ii in enumerate(np.arange(4,len(yrange),4)):
    ax = fig.add_subplot(gs[r[xx][0],r[xx][1]])
    ax.plot((xgr[ii]),(zgr[ii][0]),color='steelblue',linewidth=1)
    ax.plot((xgr[ii]),(zgr[ii][-1]),color='firebrick',linewidth=1)
    ax.grid(True)
    if (xx==0) or (xx==3) or (xx==6):
        ax.set_ylabel('Z [m NAVD88]')
    if xx>=6:
        ax.set_xlabel('X [m]')
    ax.set_xlim((-65,0))
    ax.set_ylim((-4,10))
    

ax1 = fig.add_subplot(gs[0:2, 0])
ax1.plot((xTOE_start),yrange,color='steelblue',linewidth=1,label='Pre-Storm')
ax1.plot((xTOE_end),yrange,color='firebrick',linewidth=1,label='Post-Storm')
ax1.legend(loc='best')
ax1.grid(True)

ax2 = fig.add_subplot(gs[0:2,1])
ax2.plot(cumDXTOE,yrange,color='teal',linewidth=1)
ax2.plot([0,0],[yrange[0],yrange[-1]],color='k',zorder=1)
ax2.grid(True)

ax3 = fig.add_subplot(gs[0:2,2])
ax3.plot(cumDV_net,yrange,color='gray',linewidth=1)
ax3.plot([0,0],[yrange[0],yrange[-1]],color='k',zorder=1)
ax3.grid(True)


# # # Import DRT Output File # # #
with open("DRT_20121029-20121031_36.19_-75.75.pkl",'rb') as filein:
    scenario = pickle.load(filein)

# # # Import FRF DuneLidar Data # # #
dat_class = 'geomorphology'
dat_type = 'elevationTransects'
inst_type = 'duneLidarTransect'
inst = 'duneLidarTransect'

variables = {'time', 'xFRF', 'yFRF', 'elevation', 'beachProfileQCFlag','latitude','longitude',
             'beachForeshoreSlope','lidarX','lidarY'} # geomorphology > elevationTransects > duneLidarTransect
output_dir = 'F:\PyPath\duneResponseTool'

# # # INPUT DATE RANGE DESIRED # # #
# start date for data; 
year = 2015
month = 10
day = 29
hour = 0
mint = 0
sec = 0

# end date for data
yeare = 2015
monthe = 10
daye = 31
houre = 23
minte = 59
sece = 59

# turn dates into strings
date_start = f"{year}-{month:02d}-{day:02d} {hour:02d}:{mint:02d}:{sec:02d}"
date_stop = f"{yeare}-{monthe:02d}-{daye:02d} {houre:02d}:{minte:02d}:{sece:02d}"

download_FRF_ncml(dat_class, dat_type, inst_type, date_start, date_stop, variables, output_dir)
lidar_dat = concat_ncml_daterange(date_start,date_stop,inst,variables,output_dir)

#%% # # # Compute position of MHW and dune toe (approx.) for beach slopes # # #
# cons = [0.36, 4] # NAVD88 elevation contours
# # # calc shoreline position for the lidar data
# # _, lidar_shoreline_dat = calculate_shoreline_x_position(elev,elev,xFRF,wl,wl,cons,t_lidar_dt)

# _, lidar_shoreline_dat = calculate_shoreline_x_position(lidar_dat['elevation'],lidar_dat['elevation'],lidar_dat['xFRF'],
#                                                     wl_ontolidar,wl_ontolidar,cons,pd.to_datetime(lidar_dat['time'],unit='s'))


