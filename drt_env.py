# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 09:42:49 2025

@author: sdevo
"""
import numpy as np
import pandas as pd
import netCDF4
import time as tt
import glob
import requests
import zipfile
import math
import datetime
import os 
import calendar
import csv
from scipy.interpolate import interp1d
import pytz
import tkinter as tk
from tkinter import simpledialog
import datetime

def distAway(ptY, ptX, listY, listX):
    # ptX, ptY should be a single value
    distances= np.sqrt((listX-ptX)**2 + (listY-ptY)**2)
    
    return distances

def wis_determine_node(scenario):
    # wis_determine_node: finds the closest Wave Information Studies node to
    # the provided field site latitude and longitude
    #
    # Required Inputs: 'scenario' dict variable with the following:
    #   scenario['location']['lat'] (value from -90 to 90)
    #   scenario['location']['lon'] (value from -180 to 180)

    # load metadata file that provides all usable WIS nodes for this analysis:
    wistable = pd.read_excel(f"{scenario['code_direc']}\dependencies\drt_env_station_list.xlsx",sheet_name='WIS')
    
    # find closest node to the given lat/lon
    distances = distAway(scenario['location']['lat'],scenario['location']['lon'],
                         wistable[' Lat'],wistable[' Lon'])
    minval, imin = np.nanmin(distances), np.argmin(distances)
    
    if minval > 5: # don't consider a node more than 5 degrees away
        print('No environmental node close to the selected site.')
        
    # store relevant information
    wis = {}
    wis['closest_node'] = wistable['Station'][imin]
    wis['closest_lat'] = wistable[' Lat'][imin]
    wis['closest_lon'] = wistable[' Lon'][imin]
    wis['closest_depth'] = wistable[' Depth(m)'][imin]
    wis['closest_zone'] = wistable[' Region'][imin]
    
    return wis

def websave_python(url, filename):
    """
    Downloads a file from a URL and saves it to a local file.

    Args:
        url (str): The URL of the file to download.
        filename (str): The local path and filename to save the file as.
    """
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an exception for bad status codes (4xx or 5xx)

        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"File downloaded successfully to {filename}")
    # except requests.exceptions.Timeout:
    #    print(f"Error: The request timed out after {timeout_seconds} seconds.")
        return 1
    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
        print(response.text)
        return 0
    
def webread_python(url):
        response = requests.get(url)
        try:
            decoded_content = response.content.decode('utf-8')
            cr = csv.reader(decoded_content.splitlines(), delimiter=',')
            var1 = []
            var2 = []
            for row in cr:
                var1.append(row[0])
                var2.append(row[1])
                
            var1 = var1[1:]
            var2 = var2[1:]
            
            # Convert to epoch seconds
            epoch_times = [datetime.datetime.strptime(t, '%Y-%m-%d %H:%M').replace(tzinfo=datetime.timezone.utc).timestamp() for t in var1]
            wl = []
            for t in var2: # handle spaces in data (missing data):
                t = t.strip()
                if t == '' or t.lower() == 'nan' or t.lower == 'NaN':
                    wl.append(float('nan'))
                else:
                    wl.append(float(t))
            return epoch_times, wl
        except requests.exceptions.RequestException as e:
            print(f"Error reading api data: {e}")
            print(response.text)
            
def unzipfile(zip_file_path,extraction_path):
    try:
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            # Extract all contents of the zip file to the specified directory
            zip_ref.extractall(extraction_path)
        print(f"Successfully unzipped '{zip_file_path}' to '{extraction_path}'")
    except zipfile.BadZipFile:
        print(f"Error: '{zip_file_path}' is not a valid ZIP file.")
    except FileNotFoundError:
        print(f"Error: ZIP file not found at '{zip_file_path}'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def wrapto360(angle_deg):
    """Helper for wrap_to_180 — wraps to [0, 360)."""
    angle = np.asarray(angle_deg)
    wrapped = np.mod(angle, 360)
    
    if np.isscalar(angle_deg):
        if wrapped == 0 and angle > 0:
            wrapped = 360.0
        return float(wrapped)
    else:
        wrapped[(wrapped == 0) & (angle > 0)] = 360
        return wrapped


def wrapto180(angle_deg):
    """
    Wraps angle(s) in degrees to the interval (-180, 180].
    
    Parameters
    ----------
    angle_deg : float or array-like
        Input angle(s) in degrees.
    
    Returns
    -------
    wrapped : float or ndarray
        Wrapped angle(s) in (-180, 180].
    """
    angle = np.asarray(angle_deg, dtype=float)
    wrapped = ((angle + 180) % 360) - 180

    # Ensure 180 stays 180 (not -180)
    wrapped[wrapped == -180] = 180

    # Return scalar if input was scalar
    return float(wrapped) if np.isscalar(angle_deg) else wrapped


def dispersion(w,h,showflag=True):
    # [k,n,c] = dispersion(w,h,showflag)
    # dispersion solves the linear dispersion relation given
    # frequency w (or omega) = 2pi/period [1/sec]
    # and water depth, h [m]
    # returns the wavenumber, k = 2pi/wavelength [1/m],
    # n = ratio of group speed to phase speed, and
    # c = phase speed [m/s]
    
    if h<=0 and showflag:
        print('WARNING: water depth should be > 0')
        print(('Wavenumber of 10^10 will be returned for all depths <= 0'))
        h = 1
        k = 10**10
        n = 1 
        c = 0 
    else:   
        g= 9.81
        k = w/np.sqrt(g*h)
        diff = np.max(np.max(w**2 - g*k*np.tanh(k*h)))
        while np.abs(diff)> 1*10**-8:
            knew = k - (w**2 - g*k*np.tanh(k*h))/(-g*np.tanh(k*h) - g*k*h*((1/np.cosh(k*h))**2))
            k = knew
            diff = np.max(np.max(w**2 - g*k*np.tanh(k*h)))
            
        c = w/k
        with np.errstate(over='ignore'):
            denom = np.sinh(2*k*h) # as k*h goes large, sinh --> inf (e.g., large depth)
        n = 0.5*(1+(2*k*h)/denom)
    return k, n, c
    
def snells(alpha0,h0,T,h1):
    # alpha = snells(alpha,h0,T,d)
    # Snells law to determine wave direction changes
    # calculate the changing angle due to refraction for the case
    # of straight parallel contours for monochromatic waves
    #
    # alpha0 [deg]
    # h0 [m]
    # T [s]
    # d [m] matrix or array of depths from which to find angles
    # alpha [deg] is matrix or array of angles with size(alpha)=size(d)
    omega = 2*np.pi/T
    k0 = np.full_like(omega,np.nan)
    c0 = np.full_like(omega,np.nan)
    k = np.full_like(omega,np.nan)
    c = np.full_like(omega,np.nan)
    alpha = np.full_like(omega,np.nan)
    k0, _, c0= dispersion(omega,h0)
    k, _, c = dispersion(omega,h1)
    alpha = np.real(180/np.pi * np.arcsin(np.sin(alpha0*np.pi/180)*c/c0))

    # for ww in range(len(omega)):
    #     k0[ww], _, c0[ww] = dispersion(omega[ww],h0)
    #     k[ww], _, c[ww] = dispersion(omega[ww],h1)
    #     alpha[ww] = np.real(180/np.pi * np.arcsin(np.sin(alpha0[ww]*np.pi/180)*c[ww]/c0[ww]))
    
    return alpha

def shoal_waves(H0,h0,alpha0,T,h1):
    # H1, L1,alpha1 = shoal_waves(H0, h0, alpha0, T, h1)
    #shoals surface gravity waves to a different water depth
    # H0 [m] given wave height
    # h0 [m] given water depth
    # alpha0 [deg] given angle relative to shore normal
    # T [s] wave period
    # h1 [m] depth where characteristics are sought
    
    alpha1 = snells(alpha0,h0,T,h1)
    k0 = np.full_like(T,np.nan)
    n0 = np.full_like(T,np.nan)
    c0 = np.full_like(T,np.nan)
    k1 = np.full_like(T,np.nan)
    n1 = np.full_like(T,np.nan)
    c1 = np.full_like(T,np.nan)
    k0, n0, c0 = dispersion(2*np.pi/T, h0)
    k1, n1, c1 = dispersion(2*np.pi/T, h1)
    # for ww in range(len(T)):
    #     k0[ww], n0[ww], c0[ww] = dispersion(2*np.pi/T[ww], h0)
    #     k1[ww], n1[ww], c1[ww] = dispersion(2*np.pi/T[ww], h1)
    L1 = 2*np.pi/k1
    H1 = H0*np.sqrt((c0*n0)/(c1*n1))*np.sqrt(np.cos(alpha0*np.pi/180)/np.cos(alpha1*np.pi/180))
    
    return H1, L1, alpha1

def transform_waves(waves, scenario):
    
    #fix wave problems:
    localD = np.abs(wrapto180(waves['D_deepwater']-scenario['grids']['morphometrics']['azimuth']))
    
    # also set wave heights as zero if wave direction is headed offshore
    ifind = np.where(np.abs(localD)>=90)
    Hstemp = waves['Hs_deepwater'].copy()
    Hstemp[ifind] = 0
    localD[ifind] = 0 
    
    # the refraction gets confused for really oblique waves, so set an upper limit
    maxD = 60
    ifind = np.where(np.abs(localD)>=maxD)
    Hstemp = waves['Hs_deepwater']
    localD[ifind] = maxD
    
    # shoal waves
    waves['Hs_25m'], waves['L_25m'], waves['D_25m'] = shoal_waves(np.double(Hstemp),
                                                                  np.abs(np.double(waves['depth'])),
                                                                  np.double(localD),
                                                                  np.double(waves['Tp']),
                                                                  25)
    waves['Hs_25m'] = np.real(waves['Hs_25m'])
    
    if 'Hs_25m' not in waves:
        waves['Hs_25m'] = waves['Hs_deepwater']
        waves['L_25m'] = 9.81*waves['Tp']/(2*np.pi)
        waves['D_25m'] = np.zeros(len(waves['Hs_25m']))
        
    # there are also other hiccups in the shoal code that leads to large
    # wave heights. limiter added to prevent this
    ibad = np.where(waves['Hs_25m']>20)
    waves['Hs_25m'][ibad] = waves['Hs_deepwater'][ibad]
    
    return waves
    
def wis_download(scenario):
    # wis_download: downloads WIS ONLNS file and loads data
    #
    # Required Inputs: 'scenario' dict variable with the following:
    #   scenario['wis']['closest_node']
    #   scenario['wis']['closest_zone']
    
    # relevant info from scenario file
    station_num = scenario['wis']['closest_node']
    
    # check what part of the country the zone is in
    if scenario['wis']['closest_zone'] == ' Atlantic':
        zone = 'atl' # note that there are diffferent conventions depending on whether pulling onlns or thredds
        zone2 = 'Atlantic'
    elif scenario['wis']['closest_zone'] == ' Pacific':
        zone = 'pac'
        zone2 = 'Pacific'
    elif scenario['wis']['closest_zone'] == ' GulfOfMexico':
        zone = 'gom'
        zone2 - 'GulfOfMexico'
        
    try:
        # download wave data
        url  = f"https://chlthredds.erdc.dren.mil/thredds/dodsC/wis/{zone2}/ST{station_num}/ST{station_num}.ncml#noprefetch"
        # download and convert time
        time = netCDF4.Dataset(url,'r')['time'][:]
        tt.sleep(1) # seems like need to give some time to thredds
        time = time.filled(fill_value=np.nan)
        #tunit = netCDF4.Dataset(url.'r')['time']['units'] # reading attributes of variable time
        mtime = pd.to_datetime(time,unit='s')
        
        # finding index that corresponds to dates of interest
        itime = (mtime>=pd.Timestamp(scenario['timing']['start_date'])) & (mtime<=pd.Timestamp(scenario['timing']['end_date']))
        #indices in netCDF record of data of interest
        itime = np.where(itime)[0]
        # pulling data from server with itime index
        time = mtime[itime] # record of time in datetime
        waveHs = netCDF4.Dataset(url,'r')['waveHs'][np.min(itime):np.min(itime)+len(itime)].filled(np.nan)
        #tt.sleep(1) # seems like need to give some time to thredds
        waveTp = netCDF4.Dataset(url,'r')['waveTp'][np.min(itime):np.min(itime)+len(itime)].filled(np.nan)
        #tt.sleep(1) # seems like need to give some time to thredds
        waveD = netCDF4.Dataset(url,'r')['waveMeanDirection'][np.min(itime):np.min(itime)+len(itime)].filled(np.nan)
        #tt.sleep(1) # seems like need to give some time to thredds
        windSpeed = netCDF4.Dataset(url,'r')['windSpeed'][np.min(itime):np.min(itime)+len(itime)].filled(np.nan)
        #tt.sleep(1) # seems like need to give some time to thredds
        windD = netCDF4.Dataset(url,'r')['windDirection'][np.min(itime):np.min(itime)+len(itime)].filled(np.nan)
        #tt.sleep(1) # seems like need to give some time to thredds

    except: # alternatively download text files instead using thredds server if down
        if zone2 == 'Pacific':
            #last date ==
            if scenario['timing']['end_date'] > np.datenum(2011,12,31):
                print('WIS Tredds Download Failed: Pick a Date Before 2011 for the Pacific Basin to proceed with an alternative ONLNS File')
        elif zone2 == 'Atlantic':
            if scenario['timing']['end_date'] > np.datenum(2014,12,31):
                print('WIS Tredds Download Failed: Pick a Date Before 2014 for the Atlantic Basin to proceed with an alternative ONLNS File')
        elif zone2 == 'GulfOfMexico':
            if scenario['timing']['end_date'] > np.datenum(2014,12,31):
                print('WIS Tredds Download Failed: Pick a Date Before 2011 for the Atlantic Basin to proceed with an alternative ONLNS File')
        
        try: #check to see if file is already downloaded, speeds up process substantially
            file = glob.glob(f"ST{station_num}.onlns")
            if np.size(file)==1:
                data = np.loadtxt(f"ST{station_num}.onlns")
            else:
                data = np.loadtxt(file[0])
        except: # if file not available, then download
            url_to_download = f"http://wis.usace.army.mil/data/{zone}/onlns/raw/ST{station_num}_ONLNS.zip"
            good = websave_python(url_to_download, 'wis.temp.onlns.zip')
            if not good:
                print('Error downloading wave data from WIS server.')
                
            # Specify the path to the ZIP file
            zip_file_path = 'wis.temp.onlns.zip'
            
            # Specify the directory where you want to extract the files (optional)
            # If not specified, files will be extracted to the current working directory
            extraction_path = []
            
            unzipfile(zip_file_path,extraction_path)
            file = glob.glob(f"ST{station_num}.onlns")
            if np.size(file)==1:
                data = np.loadtxt(f"ST{station_num}.onlns")
            else:
                data = np.loadtxt(file[0])

        # pull out relevant variables from the WIS ONLNS files
        time = np.datenum(data[:,0])
        _, iunique = np.unique(time) # get rid of repeat time issues here
        time = np.double(time[iunique])
        windSpeed = np.double(data[iunique,4])
        windD = np.double(data[iunique,5])
        waveHs = np.dounle(data[iunique,9])
        waveTp = np.double(data[iunique,11])
        waveD = np.double(data[iunique,15])
    
    waves = {}
    winds = {}
    # interpolate waves onto desired time interval
    waves['Hs_deepwater'] = np.interp(scenario['timing']['times'],time,waveHs)
    waves['Tp'] = np.interp(scenario['timing']['times'],time,waveTp)
    waves['D_deepwater'] = np.interp(scenario['timing']['times'],time,waveD)
    waves['depth'] = scenario['wis']['closest_depth']
    winds['windSpeed'] = np.interp(scenario['timing']['times'],time,windSpeed)
    winds['windDirection'] = np.interp(scenario['timing']['times'],time,windD)

    # shoal wave heights to 25 m for Stockdon and fix any other output issues
    waves = transform_waves(waves,scenario)
    
    return waves, winds

def noaa_determine_node(scenario):
    # noaa = noaa_determine_node(scenario)
    # finds the closest NOAA tide gauge to the provided field site lat and lon
    # Required Inputs: 'scenario' dict variable with the following:
    #   scenario['location']['lat'] (value from -90 to 90)
    #   scenario['location']['lon'] (value from -180 to 180)

    # load metadata
    noaatable = pd.read_excel(f"{scenario['code_direc']}\dependencies\drt_env_station_list.xlsx",sheet_name='NOAA_Tides')
    
    # find closest node to the given lat/lon
    distances = distAway(scenario['location']['lat'],scenario['location']['lon'],
                         noaatable['Lat'],noaatable['Lon'])
    minval, imin = np.nanmin(distances), np.argmin(distances)
    
    if minval > 5: # don't consider a node more than 5 degrees away
        print('No environmental node close to the selected site.')
        
    # store relevant information
    noaa = {}
    noaa['closest_node'] = noaatable['Station'][imin]
    noaa['closest_lat'] = noaatable['Lat'][imin]
    noaa['closest_lon'] = noaatable['Lon'][imin]
    
    return noaa

def interp1gap(*args, maxgapval=np.inf, method='linear',
               interpval=np.nan, extrap=False, extrapval=None):
    """
    Python equivalent of Chad Greene's interp1gap (MATLAB).
    Interpolates across small gaps in 1D data but leaves large gaps as NaN (or user value).
    
    Parameters
    ----------
    v : array_like
        Data vector to interpolate (if called with one argument).
    x : array_like, optional
        Independent variable. Default = index of v.
    xq : array_like, optional
        Query points. Default = same as x.
    maxgapval : float, optional
        Maximum gap over which to interpolate (in same units as x).
        Default = np.inf (interpolate over all gaps).
    method : str, optional
        Interpolation method ('linear', 'nearest', 'spline', 'cubic', 'pchip').
    interpval : float, optional
        Value to assign to points within large gaps. Default = np.nan.
    extrap : bool, optional
        Whether to allow extrapolation beyond x range.
    extrapval : float, optional
        Value to assign outside x range (if extrap=True).
        
    Returns
    -------
    vq : ndarray
        Interpolated data.
    """

    # --- Parse inputs ---
    if len(args) == 1:
        v = np.asarray(args[0])
        x = np.arange(len(v))
        xq = x.copy()
    elif len(args) == 2:
        v = np.asarray(args[0])
        maxgapval = args[1]
        x = np.arange(len(v))
        xq = x.copy()
    elif len(args) == 3:
        x, v, xq = map(np.asarray, args)
    elif len(args) == 4:
        x, v, xq, maxgapval = map(np.asarray, args)
    else:
        raise ValueError("interp1gap() accepts 1–4 positional arguments.")

    # --- Input checks ---
    if not np.isscalar(maxgapval):
        raise ValueError("maxgapval must be a scalar.")
    if v.ndim != 1 or x.ndim != 1:
        raise ValueError("x and v must be 1D vectors.")
    if xq.ndim != 1:
        raise ValueError("xq must be a 1D vector.")

    # --- Remove NaNs ---
    mask = ~np.isnan(v)
    x, v = x[mask], v[mask]

    # --- Interpolator setup ---
    if extrap:
        f = interp1d(x, v, kind=method, bounds_error=False, fill_value=extrapval)
    else:
        f = interp1d(x, v, kind=method, bounds_error=False, fill_value=np.nan)

    vq = f(xq)

    # --- Find gaps larger than maxgapval ---
    x_gap = np.diff(x)
    big_gap_indices = np.where(np.abs(x_gap) > maxgapval)[0]
    ind_int = []

    for ind in big_gap_indices:
        if x_gap[0] >= 0:  # monotonically increasing
            gap_mask = (xq > x[ind]) & (xq < x[ind + 1])
        else:  # monotonically decreasing
            gap_mask = (xq > x[ind + 1]) & (xq < x[ind])
        ind_int.extend(np.where(gap_mask)[0])

    # --- Replace values within large gaps ---
    vq[ind_int] = interpval

    return vq

def noaa_download_tides(scenario):
    #tides = noaa_download_tides(scenario)
    # function to download NOAA tides
    # Required Inputs: 'scenario' dict variable with the following:
    #   scenario['timing']['start_date'] (model start date in datenum format)
    #   scenario['timing']['end_date'] (model end date in datenum format)
    #   scenario['noaa']['closest_node'] 
    
    # initialize variables
    startYear = scenario['timing']['start_date'].year
    endYear = scenario['timing']['end_date'].year
    datum = 'NAVD' # stations in the list should all be NAVD compatible
    gauge = scenario['noaa']['closest_node']
    
# loop through each year of verified tide data for noaa erddap server
    wl = []
    time = []
    for yr in np.arange(startYear,endYear+1,1):
        #erdap server link:
        website = f'https://opendap.co-ops.nos.noaa.gov/erddap/tabledap/IOOS_Hourly_Height_Verified_Water_Level.nc?STATION_ID%2CDATUM%2CBEGIN_DATE%2CEND_DATE%2Ctime%2CWL_VALUE&STATION_ID=%22{gauge}%22&DATUM=%22{datum}%22&BEGIN_DATE%3E=%22{yr}0101%2000%3A00%22&END_DATE%3C=%22{yr}0131%2023%3A59%22'
        output_name = 'tides.nc'
        # options = {}
        # options['timeout'] = 120
        good = websave_python(website, output_name)
        if good:
            # load and then clean up variables
            file2read = netCDF4.Dataset(output_name,'r')
            
            #store data
            wltemp = file2read.variables['WL_VALUE'][:].filled(np.nan)
            timetemp = file2read.variables['time'][:].filled(np.nan)
            # timetemp = double(timetemp/86400 + datenum(1970,1,1));
            wl.append(wltemp)#[wl; wltemp]
            time.append(timetemp)#[time; timetemp];
            file2read.close()
            os.remove(output_name)
        else:
            print('Error downloading from NOAA ERDDAP. Moving to NOAA API Download.')
            #read directly from api (issues sometimes with NAVD datum from erddap server?)
            for im in np.arange(1,13):
                st_date = datetime.date(yr, im, 1)
                # End date = last day of the month
                last_day = calendar.monthrange(yr, im)[1]
                end_date = datetime.date(yr, im, last_day)
                st_date_str = st_date.strftime('%Y%m%d')
                end_date_str = end_date.strftime('%Y%m%d')
                url= f'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=water_level&application=NOS.COOPS.TAC.WL&begin_date={st_date_str}&end_date={end_date_str}&datum={datum}&station={gauge}&time_zone=GMT&units=metric&interval=h&format=CSV'
                
                timetemp, wltemp = webread_python(url)
                
                wl.append(wltemp)#[wl; wltemp]
                time.append(timetemp)#[time; timetemp];
    wl = np.concatenate(wl)
    time = np.concatenate(time)
    del wltemp, timetemp
        
    # loop through each year of predicted tide data to fill in any data gaps in ther verified data
    yrs = np.arange(startYear,endYear+1) #define all the year withing the start and end year limit
    wl_pred = []
    time_pred = []
    for i in np.arange(0,len(yrs)):
        st_date = datetime.datetime(yrs[i],1,1)
        end_date = datetime.datetime(yrs[i],12,31)
        st_date_str = st_date.strftime('%Y%m%d')
        end_date_str = end_date.strftime('%Y%m%d')
        url=f'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date={st_date_str}&end_date={end_date_str}&datum={datum}&station={gauge}&time_zone=GMT&units=metric&interval=h&format=CSV'
       
        timetemp, wltemp = webread_python(url)

        wl_pred.append(wltemp)#[wl; wltemp]
        time_pred.append(timetemp)#[time; timetemp];
    
    wl_pred = np.concatenate(wl_pred)
    time_pred = np.concatenate(time_pred)
    del wltemp, timetemp
    
    scenario_epoch = (scenario['timing']['times'] - pd.Timestamp("1970-01-01"))/pd.Timedelta("1s")
    # interpolate data to output
    tides_wl = interp1gap(time,wl,scenario_epoch, maxgapval=6*3600) # max gap is 6 hours  
    f = interp1d(time_pred, wl_pred, kind='linear', bounds_error=False, fill_value=np.nan)
    tides_pwl = f(scenario_epoch)  
    ibad = np.where(np.isnan(tides_wl))[0]
    tides_wl[ibad] = tides_pwl[ibad]   

    tides = {}
    tides['wl'] = tides_wl
    tides['pwl'] = tides_pwl
    return tides
  
def noaa_download_tides_prediction(scenario):

    #initialize variables
    startYear = scenario['timing']['start_date'].year
    endYear = scenario['timing']['end_date'].year
    datum = 'NAVD' # stations in the list should all be NAVD compatible
    gauge = scenario['noaa']['closest_node']

    yrs = np.arange(startYear,endYear+1) #define all the year withing the start and end year limit
    wl_pred = []
    time_pred = []
    for i in np.arange(0,len(yrs)):
        st_date = datetime.datetime(yrs[i],1,1)
        end_date = datetime.datetime(yrs[i],12,31)
        st_date_str = st_date.strftime('%Y%m%d')
        end_date_str = end_date.strftime('%Y%m%d')
        url=f'https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=predictions&application=NOS.COOPS.TAC.WL&begin_date={st_date_str}&end_date={end_date_str}&datum={datum}&station={gauge}&time_zone=GMT&units=metric&interval=h&format=CSV'
       
        timetemp, wltemp = webread_python(url)

        wl_pred.append(wltemp)#[wl; wltemp]
        time_pred.append(timetemp)#[time; timetemp];
    
    wl_pred = np.concatenate(wl_pred)
    time_pred = np.concatenate(time_pred)
    del wltemp, timetemp

    try:
        scenario_epoch = (scenario['timing']['times'] - pd.Timestamp("1970-01-01"))/pd.Timedelta("1s")
    except:
        scenario_epoch = scenario['timing']['times']    
    # interpolate data to output
    f = interp1d(time_pred, wl_pred, kind='linear', bounds_error=False, fill_value=np.nan)
    tides_pwl = f(scenario_epoch)  
    tides = {}
    tides['wl'] = tides_pwl
    return tides

def download_ww3(start_time, end_time, lat, lon):
    # waves = download_ww3(start_time, end_time, lat, lon)
    #download waves from WaveWatchIII
    try:
        start_str = start_time.strftime("%Y-%m-%dT00:00:00Z")
        end_str = end_time.strftime("%Y-%m-%dT00:00:00Z")
    except:
        start_str = start_time
        end_str = end_time
    lonn = wrapto360(lon)
    url = f"https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ww3_global.nc?Tdir%5B({start_str}):1:({end_str})%5D%5B(0.0):1:(0.0)%5D%5B({lat}):1:({lat})%5D%5B({lonn}):1:({lonn})%5D,Tper%5B({start_str}):1:({end_str})%5D%5B(0.0):1:(0.0)%5D%5B({lat}):1:({lat})%5D%5B({lonn}):1:({lonn})%5D,Thgt%5B({start_str}):1:({end_str})%5D%5B(0.0):1:(0.0)%5D%5B({lat}):1:({lat})%5D%5B({lonn}):1:({lonn})%5D"
    # url = "https://pae-paha.pacioos.hawaii.edu/erddap/griddap/ww3_global.nc?Tdir%5B(2025-10-22T00:00:00Z):1:(2025-10-26T00:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(36.25):1:(36.25)%5D%5B(286.42):1:(286.42)%5D,Tper%5B(2025-10-22T00:00:00Z):1:(2025-10-26T00:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(36.25):1:(36.25)%5D%5B(286.42):1:(286.42)%5D,Thgt%5B(2025-10-22T00:00:00Z):1:(2025-10-26T00:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(36.25):1:(36.25)%5D%5B(286.42):1:(286.42)%5D"
    # # url = f"https://coastwatch.pfeg.noaa.gov/erddap/griddap/NWW3_Global_Best.csv?Tdir%5B({start_str}):1:({end_str})%5D%5B(0.0):1:(0.0)%5D%5B({lat}):1:({lat})%5D%5B({wrapto360(lon)}):1:({wrapto360(lon)})%5D,Tper%5B({start_str}):1:({end_str})%5D%5B(0.0):1:(0.0)%5D%5B({lat}):1:({lat})%5D%5B({wrapto360(lon)}):1:({wrapto360(lon)})%5D,Thgt%5B({start_str}):1:({end_str})%5D%5B(0.0):1:(0.0)%5D%5B({lat}):1:({lat})%5D%5B({wrapto360(lon)}):1:({wrapto360(lon)})%5D"
    good = websave_python(url, 'ww3.nc')    
    
    if not good:
        print('Error downloading WW3 Forecast.')
        # artificial nan fill:
        waves = {}
        waves['Hs_deepwater'] = np.full(10,np.nan)
        waves['time'] = np.full(10,np.nan)
        waves['D_deepwater'] = np.full(10,np.nan)
        waves['Tp'] = np.full(10,np.nan)
        waves['latitude'] = np.full(10,np.nan)
        waves['longitude'] = np.full(10,np.nan)
    else:
        # load and then clean up variables
        file2read = netCDF4.Dataset('ww3.nc','r')
        
        #store data
        waves = {}
        waves['Hs_deepwater'] = file2read.variables['Thgt'][:,0,0,0].filled(np.nan)
        waves['time'] = file2read.variables['time'][:].filled(np.nan)
        waves['D_deepwater'] = file2read.variables['Tdir'][:,0,0,0].filled(np.nan)
        waves['Tp'] = file2read.variables['Tper'][:,0,0,0].filled(np.nan)
        waves['latitude'] = file2read.variables['latitude'][:].filled(np.nan)
        waves['longitude'] = file2read.variables['longitude'][:].filled(np.nan)
        file2read.close()
        os.remove('ww3.nc')
    return waves

def ww3_forecast_download(scenario):
    #waves2, tides = ww3_forecast_download(scenario):
    # downloads wave watch III forecast
    # Required Inputs: 'scenario' dict variable with the following:
    #       scenario.timing.start_date
    #       scenario.timing.end_date     
    #       scenario.location.lat [value from -90 to 90]
    #       scenario.location.lon [value from -180 to 180]       

    # generate timing in format needed for erddap
    start_time = (datetime.datetime.utcnow().replace(hour=0, minute=0, second=0, microsecond=0))#.strftime("%Y-%m-%dT00:00:00Z")
    end_time = (datetime.datetime.utcnow().replace(hour=0, minute=0, second=0, microsecond=0) 
                + datetime.timedelta(days=1))#.strftime("%Y-%m-%dT00:00:00Z")
    # first just try the coordinates given
    waves = download_ww3(start_time, end_time, scenario['location']['lat'], scenario['location']['lon'])
    
    # next try the location of the closest wis node
    if len(np.where(np.isnan(waves['Hs_deepwater'][:]))[0]) > 0:
        wis = wis_determine_node(scenario)
        lat= wis['closest_lat']
        lon = wis['closest_lon']
        waves = download_ww3(start_time, end_time, lat,lon);
   
    # now try different locations if not
    # go west if in the west coast, east on the east coast, and south in the
    # gulf coast
    if len(np.where(np.isnan(waves['Hs_deepwater'][:]))[0]) > 0:
        if wis['closest_zone'] == ' Atlantic':
            lon+=1
        elif wis['closest_zone'] == ' Pacific':
            lon-=1;
        elif wis['closest_zone'] == 'GulfOfMexico':
            lat-=1    
        waves = download_ww3(start_time, end_time, lat,lon) 

    # try one last time
    if len(np.where(np.isnan(waves['Hs_deepwater'][:]))[0]) > 0:
        if wis['closest_zone'] == ' Atlantic':
            lon+=1
        elif wis['closest_zone'] == ' Pacific':
            lon-=1;
        elif wis['closest_zone'] == 'GulfOfMexico':
            lat-=1    
        waves = download_ww3(start_time, end_time, lat,lon) 
   
    # download local water depths at wave note since this is not part of the WW3 outout
    url = f"https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm15plus.nc?z%5B({waves['latitude'][0]}):1:({waves['latitude'][0]})%5D%5B({waves['longitude'][0]-360}):1:({waves['longitude'][0]-360})%5D"
            # https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm15plus.nc?z%5B(-90.0):1:(90.0)%5D%5B(-180.0):1:(180.0)%5D
    good = websave_python(url,'srtm.nc')
    
    #load water depth data
    file2read = netCDF4.Dataset('srtm.nc','r')
    waves['depth'] = file2read.variables['z'][:,0].filled(np.nan)
    file2read.close()
    os.remove('srtm.nc')

    #transform waves to shallow water
    waves = transform_waves(waves, scenario)
    
    #generate new timings
    scenario['timing']['times'] = np.arange(np.round(np.min(waves['time'])),np.floor(np.max(waves['time'])),scenario['timing']['dt']*3600) # if dt is in hours, since waves['time'] in seconds (epoch time)
   
    # interpolate
    waves2 = {}
    waves2['Hs_25m'] = np.interp(scenario['timing']['times'],waves['time'], waves['Hs_25m'])
    waves2['Hs_deepwater'] = np.interp(scenario['timing']['times'],waves['time'],waves['Hs_deepwater'])
    waves2['L_25m'] = np.interp(scenario['timing']['times'],waves['time'], waves['L_25m'])
    waves2['D_25m'] = np.interp(scenario['timing']['times'],waves['time'], waves['D_25m'])
    waves2['D_deepwater'] = np.interp(scenario['timing']['times'],waves['time'], waves['D_deepwater'])
    waves2['Tp'] = np.interp(scenario['timing']['times'],waves['time'], waves['Tp'])
    waves2['times'] = scenario['timing']['times']  
    times = scenario['timing']['times']  
    return  waves2, times

def uv_to_wswd(u, v):
    # ws, wd = uv_to_wswd(u,v)
    #function that takes vectors (u,v) and converts to  
    # wind speed and wind direction

    ws = np.full(len(u),np.nan)
    wd = np.full(len(u),np.nan)
    e = np.where((~np.isnan(u)) & (~np.isnan(v)))[0]
    ws[e] = np.sqrt(u[e]*u[e] + v[e]*v[e])
    wd[e] = 270 - (180/np.pi)*np.arctan2(v[e],u[e])
    
    for i in range(len(wd)):
        if wd[i] > 360:
            wd[i]-=360
    return ws, wd
 
def gfs_forecast_download(scenario):
    # winds = gfs_forecast_download(scenario)
    
    # generate timing in format needed for erddap
    start_time = (datetime.datetime.utcnow().replace(hour=0, minute=0, second=0, microsecond=0))#.strftime("%Y-%m-%dT00:00:00Z")
    end_time = (datetime.datetime.utcnow().replace(hour=0, minute=0, second=0, microsecond=0) 
                + datetime.timedelta(days=1))#.strftime("%Y-%m-%dT00:00:00Z")

    try:
        start_str = start_time.strftime("%Y-%m-%dT00:00:00Z")
        end_str = end_time.strftime("%Y-%m-%dT00:00:00Z")
    except:
        start_str = start_time
        end_str = end_time
    lonn = wrapto360(scenario['location']['lon'])
    url = f"https://coastwatch.pfeg.noaa.gov/erddap/griddap/NCEP_Global_Best.nc?ugrd10m%5B({start_str}):1:({end_str})%5D%5B({scenario['location']['lat']}):1:({scenario['location']['lat']})%5D%5B({lonn}):1:({lonn})%5D,vgrd10m%5B({start_str}):1:({end_str})%5D%5B({scenario['location']['lat']}):1:({scenario['location']['lat']})%5D%5B({lonn}):1:({lonn})%5D"
    good = websave_python(url,'gfs.nc')
    
    #load wind forecast data
    file2read = netCDF4.Dataset('gfs.nc','r')
    time = file2read.variables['time'][:].filled(np.nan)
    u = file2read.variables['ugrd10m'][:,0,0].filled(np.nan)
    v = file2read.variables['vgrd10m'][:,0,0].filled(np.nan)
    windSpeed, windD = uv_to_wswd(u, v)
   
    # interpolate waves onto desired time interval
    winds = {}
    winds['windSpeed'] = np.interp(scenario['timing']['times'],time, windSpeed)
    winds['windDirection'] = np.interp(scenario['timing']['times'],time, windD)
    file2read.close()
    os.remove('gfs.nc')
    
    return winds

def download_ESTOFS(scenario, zone):
    # tides = download_ESTOFS(scenario,zone)
    #Function to download ESTOFS surge data for input to the model
    input_lat = scenario['location']['lat']
    input_lon = scenario['location']['lon']
            
    # Time keeping
   # Current time
    time = datetime.datetime.now()
    hour = time.hour
    
    # If the current hour >= 21 (i.e., 9 PM), use tomorrow's date
    if hour >= 21:
        datestring = (time + datetime.timedelta(days=1)).strftime('%Y%m%d')
    else:
        datestring = time.strftime('%Y%m%d')

    # Determine which geographic zone are in
    if zone == ' Pacific':
        estofs_zone = 'west'
    elif zone == ' Atlantic':
        estofs_zone = 'east'
    elif zone == 'GulfOfMexico':
        estofs_zone = 'east'
    else:
        print('No forecast data in this zone')

    # # # Use only in forecast mode
    # # # file = ['https://nomads.ncep.noaa.gov:9090/dods/estofs_', estofs_zone,'/',datestring,'/estofs_', estofs_zone,'_conus_00z'];
    # # file = f"https://nomads.ncep.noaa.gov/dods/estofs/{datestring}/estofs_conus.{estofs_zone}_00z"
    # # try:
    # #     dataset = netCDF4.Dataset(file)
    # #     print("Variables available:", list(dataset.variables.keys()))
    # #     # dataset.close()
    # # except:        
    # #     try:
    # #         datestring = datetime.datetime.now().strftime('%Y%m%d')
    # #         file = f"https://nomads.ncep.noaa.gov/dods/estofs/{datestring}/estofs_conus.{estofs_zone}_00z"
    # #         dataset = netCDF4.Dataset(file)
    # #         print("Variables available:", list(dataset.variables.keys()))
    # #         # dataset.close()
    # #     except:
    # #         datestring = (datetime.datetime.now() - datetime.timedelta(days=1)).strftime('%Y%m%d')
    # #         file = f"https://nomads.ncep.noaa.gov/dods/estofs/{datestring}/estofs_conus.{estofs_zone}_00z"
    # #         dataset = netCDF4.Dataset(file)
    # #         print("Variables available:", list(dataset.variables.keys()))
    # #         # dataset.close()
               
    # # Use only in forecast mode
    # # file = ['https://nomads.ncep.noaa.gov:9090/dods/estofs_', estofs_zone,'/',datestring,'/estofs_', estofs_zone,'_conus_00z'];
    # file = f"https://nomads.ncep.noaa.gov/dods/stofs_2d_glo/{datestring}/stofs_2d_glo_conus.{estofs_zone}_00z"
    # try:
    #     dataset = netCDF4.Dataset(file)
    #     print("Variables available:", list(dataset.variables.keys()))
    #     # dataset.close()
    # except:        
    #     try:
    #         datestring = datetime.datetime.now().strftime('%Y%m%d')
    #         file = f"https://nomads.ncep.noaa.gov/dods/stofs_2d_glo/{datestring}/stofs_2d_glo_conus.{estofs_zone}_00z"
    #         dataset = netCDF4.Dataset(file)
    #         print("Variables available:", list(dataset.variables.keys()))
    #         # dataset.close()
    #     except:
    #         datestring = (datetime.datetime.now() - datetime.timedelta(days=1)).strftime('%Y%m%d')
    #         file = f"https://nomads.ncep.noaa.gov/dods/stofs_2d_glo/{datestring}/stofs_2d_glo_conus.{estofs_zone}_00z"
    #         dataset = netCDF4.Dataset(file)
    #         print("Variables available:", list(dataset.variables.keys()))
    #         # dataset.close()
        
            
    # file = f"https://nomads.ncep.noaa.gov:9090/dods/estofs_{estofs_zone}/{datestring}/estofs_{estofs_zone}_conus_00z"
    london = pytz.timezone("Europe/London")
    now_london = datetime.datetime.now(london)
    
    # Equivalent to datetime('now','TimeZone','Europe/London') - 4/24
    adjusted_time = now_london - datetime.timedelta(hours=4)
    datestring = (adjusted_time - datetime.timedelta(hours=adjusted_time.hour % 24,
                                            minutes=adjusted_time.minute,
                                            seconds=adjusted_time.second)).strftime("%Y%m%d")
    # file = f"https://nomads.ncep.noaa.gov/dods/estofs/{datestring}/estofs_conus.{estofs_zone}_18z"
    # try:
    #     dataset = netCDF4.Dataset(file)
    #     print("Variables available:", list(dataset.variables.keys()))
    #     # dataset.close()
    # except:
    #     try:
    #         datestring = datetime.datetime.now().strftime('%Y%m%d')
    #         file = f"https://nomads.ncep.noaa.gov/dods/estofs/{datestring}/estofs_conus.{estofs_zone}_12z"
    #         dataset = netCDF4.Dataset(file)
    #         print("Variables available:", list(dataset.variables.keys()))
    #         # dataset.close()        
    #     except:
    #         try:
    #             datestring = (datetime.datetime.now() - datetime.timedelta(days=1)).strftime('%Y%m%d')
    #             file = f"https://nomads.ncep.noaa.gov/dods/estofs/{datestring}/estofs_conus.{estofs_zone}_06z"
    #             dataset = netCDF4.Dataset(file)
    #             print("Variables available:", list(dataset.variables.keys()))
    #             # dataset.close()  
    #         except:
    #             datestring = (datetime.datetime.now() - datetime.timedelta(days=1)).strftime('%Y%m%d')
    #             file = f"'https://nomads.ncep.noaa.gov/dods/estofs/{datestring}/estofs_conus.{estofs_zone}_00z"
    #             dataset = netCDF4.Dataset(file)
    #             print("Variables available:", list(dataset.variables.keys()))
    #             # dataset.close() 
    file = f"https://nomads.ncep.noaa.gov/dods/stofs_2d_glo/{datestring}/stofs_2d_glo_conus.{estofs_zone}_18z"
    try:
        dataset = netCDF4.Dataset(file)
        print("Variables available:", list(dataset.variables.keys()))
        # dataset.close()
    except:
        try:
            datestring = datetime.datetime.now().strftime('%Y%m%d')
            file = f"https://nomads.ncep.noaa.gov/dods/stofs_2d_glo/{datestring}/stofs_2d_glo_conus.{estofs_zone}_12z"
            dataset = netCDF4.Dataset(file)
            print("Variables available:", list(dataset.variables.keys()))
            # dataset.close()        
        except:
            try:
                datestring = (datetime.datetime.now() - datetime.timedelta(days=1)).strftime('%Y%m%d')
                file = f"https://nomads.ncep.noaa.gov/dods/stofs_2d_glo/{datestring}/stofs_2d_glo_conus.{estofs_zone}_06z"
                dataset = netCDF4.Dataset(file)
                print("Variables available:", list(dataset.variables.keys()))
                # dataset.close()  
            except:
                datestring = (datetime.datetime.now() - datetime.timedelta(days=1)).strftime('%Y%m%d')
                file = f"'https://nomads.ncep.noaa.gov/dods/stofs_2d_glo/{datestring}/stofs_2d_glo_conus.{estofs_zone}_00z"
                dataset = netCDF4.Dataset(file)
                print("Variables available:", list(dataset.variables.keys()))
                # dataset.close() 
                            
    # find the closest node
    attrs = dataset.variables['lon'].__dict__  # or adjust variable name/index
    min = attrs['minimum']      
    step = attrs['resolution']
    max = attrs['maximum']
    lon = np.double(np.arange(min,max,step))
    attrs = dataset.variables['lat'].__dict__  # or adjust variable name/index
    min = attrs['minimum']      
    step = attrs['resolution']
    max = attrs['maximum']
    lat = np.double(np.arange(min,max,step))
    LON, LAT = np.meshgrid(lon,lat)
    
    minval, ilatuse = np.min(np.abs(lat - input_lat)), np.argmin(np.abs(lat - input_lat))
    minval, ilonuse = np.min(np.abs(lon - input_lon)), np.argmin(np.abs(lon - input_lon))
    
    # set up timing information
    attrs = dataset.variables['time'].__dict__  # or adjust variable name/index
    min = attrs['minimum']      
    step = attrs['resolution']
    max = attrs['maximum']
    temptime1 = min
    hour = int(temptime1[:2])  # hour
    date_str = temptime1[3:]   # 'ddmmmyyyy'
    dt = datetime.datetime.strptime(date_str, "%d%b%Y")
    temptime1 = dt + datetime.timedelta(hours=hour)
        
    temptime2 = max
    hour = int(temptime2[:2])  # hour
    date_str = temptime2[3:]   # 'ddmmmyyyy'
    dt = datetime.datetime.strptime(date_str, "%d%b%Y")
    temptime2 = dt + datetime.timedelta(hours=hour)    
    
    dt = step
    time = pd.date_range(temptime1,temptime2,freq='1h').astype('int64') // 10**9
    
    if (dt> 0.04) and (dt<0.045):
        dt = 1/24
    
    # find longitude location closest to shore (but must be at or west of defined point) with surge values
    surge2 = np.nan
    count = 0
    while np.sum(np.isnan(surge2))>0:
        #surge = ncread(['https://nomads.ncep.noaa.gov:9090/dods/estofs_', estofs_zone, '/',datestring,'/estofs_', estofs_zone,'_conus_00z'], 'etsrgsfc', [ilonuse+count ilatuse 1], [1 1 Inf]);
        surge = dataset.variables['etsrgsfc'][:, ilatuse, ilonuse+count].filled(np.nan)
        del surge2
        surge2 = surge[:]
        count-=1
        if count<-10:
            surge2 = np.zeros(len(surge2))
            
    tides = {}
    tides['surge'] = np.interp(scenario['timing']['times'],time, surge2)
    tides['wl'] = scenario['env']['tides']['wl'] + tides['surge']
   
    return tides

def choosedialog_surge():
    options = ["No SLA", "0.10 m SLA", "0.25 m SLA", "0.50 m SLA", "1.0 m SLA"]
    values = {"No SLA":0, "0.10 m SLA":0.1, "0.25 m SLA":0.25, "0.50 m SLA":0.5, "1.0 m SLA":1.0}

    surge = [0]  # use a mutable object to store result

    def submit():
        surge[0] = values[var.get()]
        root.destroy()

    root = tk.Tk()
    root.title("Select One")
    tk.Label(root, text="Surge Forecast Download Failed.\nManually add a sea level anomaly?").pack(pady=10)

    var = tk.StringVar(root)
    var.set(options[0])  # default values
    tk.OptionMenu(root, var, *options).pack(pady=5)

    tk.Button(root, text="Close", command=submit).pack(pady=10)

    root.mainloop()
    return surge[0]

def drt_search_morphology(scenario, morph_file='drt_morphology.xlsx', sheet_name=0):
    """
    Find the nearest morphology data point to the scenario location
    and assign its dune parameters into scenario['grids']['morphometrics'].

    """

    # Load the morphology dataset
    all_morph = pd.read_excel(morph_file, sheet_name=sheet_name)

    lat_target = scenario['location']['lat']
    lon_target = scenario['location']['lon']

    # Compute Euclidean distance in degrees
    distances = np.sqrt((all_morph[' Lon'] - lon_target)**2 + (all_morph[' Lat'] - lat_target)**2)

    # Find index of closest point
    imin = np.nanargmin(distances)
    minval = distances.iloc[imin]

    # Check distance threshold (degrees)
    if minval > 1:
        print("No pre-loaded morphology data close to the selected site.")
        return scenario  # stop if too far away

    # Initialize the grids dict if not present
    scenario.setdefault('grids', {})
    scenario['grids'].setdefault('morphometrics', {})


    # Assign morphology values from the closest record
    scenario['grids']['morphometrics']['dhigh'] = float(all_morph.loc[imin, 'DuneCrestElev_m_navd'])
    scenario['grids']['morphometrics']['dtoe'] = float(all_morph.loc[imin, 'DuneToeElev_m_navd'])
    scenario['grids']['morphometrics']['duneslope'] = float(all_morph.loc[imin, 'DuneSlope'])
    scenario['grids']['morphometrics']['backshoreslope'] = float(all_morph.loc[imin, 'BeachSlope'])
    scenario['grids']['morphometrics']['azimuth'] = float(all_morph.loc[imin, 'Azimuth'])

    # Cap duneslope at 1 (angle of repose limit)
    if scenario['grids']['morphometrics']['duneslope'] > 1:
        scenario['grids']['morphometrics']['duneslope'] = 1.0

    print(f"Closest morphology record found to input Lat = {lat_target} Lon = {lon_target}: "
          f"Lat={all_morph.loc[imin, ' Lat']:.3f}, Lon={all_morph.loc[imin, ' Lon']:.3f}")
    print(f"   Dune crest elev: {scenario['grids']['morphometrics']['dhigh']} m NAVD")
    print(f"   Dune toe elev: {scenario['grids']['morphometrics']['dtoe']} m NAVD")
    print(f"   Dune slope: {scenario['grids']['morphometrics']['duneslope']} | Azimuth: {scenario['grids']['morphometrics']['azimuth']}°")

    return scenario

def drt_env(scenario):
    # drt_env: code to download relevant environmental forcings and set
    # these time series datasets needed to run erosion/accretion model
    #
    # Required Inputs: 'scenario' dict variable with the following:
    #   scenario['location']['lat'] (value from -90 to 90)
    #   scenario['location']['lon'] (value from -180 to 180)
    #   scenario['timing']['start_date'] (model start date in datenum format)
    #   scenario['timing']['end_date'] (model end date in datenum format)
    #   scenario['timing']['dt'] (model time step in hrs)
    #   scenario['type'] (string variable indicating this is a 'hindcast' or 'forecast')
    #
    # Outputs:
    #   scenario['env']['waves']
    #   scenario['env']['winds']
    #   scenario['env']['tides']
    
    if scenario['type'] == 'hindcast': #workflow for hindcast model setup
    # throw an error if the end date selected is past the period where WIS waves are present
        if scenario['timing']['end_date'] > datetime.datetime(2023,12,31):
            print('Model simulation dates must be between 1980 and 2023. Select an end date thht falls in this range')
        
        # develop scenario timing
        scenario['timing']['times'] = pd.date_range(start=scenario['timing']['start_date'],
                                                    end=scenario['timing']['end_date'],
                                                    freq=f"{scenario['timing']['dt']/24}D")

        # acquire wind and wave data from the Wave Information Studies database
        scenario['wis'] = wis_determine_node(scenario) # find closest WIS node to the relevant site
        
        # download data
        scenario['env']['waves'], scenario['env']['winds'] = wis_download(scenario)
        
        # acquire tide data for NOAA Tides and Currents
        scenario['noaa'] = noaa_determine_node(scenario)
        # scenario['env']['tides'] = noaa_download_tides(scenario)
        
        try: #download the verified tide data, if avaialble
            scenario['env']['tides'] = noaa_download_tides(scenario)
        except: #add in an error check an error can occur if no measured data exists during the time period wanted
            scenario['env']['tides'] = noaa_download_tides_prediction(scenario)
            
    elif scenario['type'] == 'forecast':  # workflow for forecast model setup

        # first pull the wave data from the global WaveWatchIII model
        scenario['env']['waves'], scenario['timing']['times'] = ww3_forecast_download(scenario)

        # create timing information based on output from WaveWatchIII
        scenario['timing']['start_date'] = pd.to_datetime(np.min(scenario['timing']['times']),unit='s').date()
        scenario['timing']['end_date'] = pd.to_datetime(np.max(scenario['timing']['times']),unit='s').date()
        
        # acquire wind data from the Global Forecast System
        scenario['env']['winds'] = gfs_forecast_download(scenario)
       
        # download tidal data from NOAA Tides and Currents
        scenario['noaa'] = noaa_determine_node(scenario)
        scenario['env']['tides'] = noaa_download_tides_prediction(scenario)

        # add in storm surge
        scenario['wis'] = wis_determine_node(scenario) # will not use wis data, but this easily tells us which coast the site is on
        zone = scenario['wis']['closest_zone']
        
        try: # the OpenDap ESTOFs site seems to have dependability issues, so a try catch is necessary here
            scenario['env']['tides'] = download_ESTOFS(scenario, zone)
        except:
            scenario['env']['tides']['surge'] = choosedialog_surge()
            scenario['env']['tides']['wl'] = scenario['env']['tides']['wl'] + scenario['env']['tides']['surge']

    return scenario