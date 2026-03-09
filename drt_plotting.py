# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 11:03:28 2025

@author: sdevo
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import gridspec

def drt_plotting(scenario):
    # drt_plotting: code to plot basic model output for deterministic run  
    # fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    # fig.subplots_adjust(wspace=0.25, hspace=0.35)
    fig = plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(2, 3, figure=fig, wspace=0.25, hspace=0.35)
    
    # Top row (3 subplots)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    
    # Bottom row:
    ax5 = fig.add_subplot(gs[1, 0:2])  # spans bottom-left two columns
    ax4 = fig.add_subplot(gs[1, 2])    # bottom-right panel
    # figure('units','normalized','outerposition',[0 0 1 1])
    
    # ---- COMMON time limits ----
    times = scenario['timing']['times']
    xlims = [np.min(times), np.max(times)]
    xpts = np.linspace(xlims[0].toordinal(), xlims[1].toordinal(), 5)
    xpts = np.unique(np.floor(xpts))
    
    # Format x-axis for date ticks
    date_fmt = mdates.DateFormatter('%Y-%m-%d')

    # # Subplot for Wind Time Series
    # ax1 = axs[0,0]
    ax1.plot(times, scenario['env']['winds']['windSpeed'],'k-',linewidth=3)
    ax1.set_ylabel('Wind Speed [m/s]')
    ax1.grid(True)
    ax1.set_title('Winds', fontweight='bold')
    ax1.set_xlim(xlims)
    ax1.xaxis.set_major_formatter(date_fmt)
    ax1.tick_params(axis='x',rotation=30)
    ax1.tick_params(width=1.5)
    ax1.spines[:].set_linewidth(1.5)
    
    ax1b = ax1.twinx()
    ax1b.plot(times,scenario['env']['winds']['windDirection'],'.',linewidth=0.5,color=[0.5, 0.5, 0.5],markersize=6) 
    ax1b.set_ylabel(r'Wind Direction $[^{\circ}]$')
    ax1b.set_ylim([0,360])
    ax1.yaxis.label.set_color('k')
    ax1b.yaxis.label.set_color([0.5, 0.5, 0.5])
    
    # # Subplot for Wave Time Series   
    # ax2 = axs[0,1]
    ax2.plot(times, scenario['env']['waves']['Hs_deepwater'],'k-',linewidth=3)
    ax2.set_ylabel(r'$H_s$ [m]')
    ax2.grid(True)
    ax2.set_title('Waves',fontweight='bold')
    ax2.set_xlim(xlims)
    ax2.xaxis.set_major_formatter(date_fmt)
    ax2.tick_params(axis='x',rotation=30)
    ax2.spines[:].set_linewidth(1.5)
    
    ax2b  = ax2.twinx()
    ax2b.plot(times,scenario['env']['waves']['D_deepwater'],'.',color=[0.5, 0.5, 0.5],markersize=6)
    ax2b.set_ylabel(r'Wave Direction $[^{\circ}]$')
    ax2b.set_ylim([0,360])
    ax2.yaxis.label.set_color('k')
    ax2b.yaxis.label.set_color([0.5, 0.5, 0.5])
    
    # # Subplot for Water Level Time Series
    # ax3 = axs[0,2]
    ax3.plot(times, scenario['env']['tides']['wl'],color=[0.5, 0.5, 0.5], linewidth=3,label='SWL')
    ax3.plot(times,scenario['erosion']['TWL'],'k-',linewidth=3,label='TWL')
    
    dtoe = scenario['grids']['morphometrics']['dtoe']
    dhigh = scenario['grids']['morphometrics']['dhigh']
    ax3.plot(xlims,[dtoe, dtoe],'--',color=[0.9, 0.3, 0.3],linewidth=2)
    ax3.plot(xlims,[dhigh, dhigh],'--',color=[0.3, 0.1, 0.1],linewidth=2)
    xo = xlims[0] + (xlims[1]-xlims[0])*0.05
    ax3.text(xo, dhigh+0.25,'Dune Crest',fontweight='bold',color=[0.3, 0.1, 0.1])
    ax3.text(xo, dtoe+0.25,'Dune Toe',fontweight='bold',color=[0.9, 0.3, 0.3])
    
    ax3.set_ylabel('SWL & TWL [m]')
    ax3.set_title('Water Levels',fontweight='bold')
    ax3.grid(True)
    ax3.set_xlim(xlims)
    ax3.spines[:].set_linewidth(1.5)
    ax3.xaxis.set_major_formatter(date_fmt)
    ax3.tick_params(axis='x',rotation=30)
    ylims = [np.min((np.min(scenario['env']['tides']['wl']),dtoe)) - 1, np.max((np.max(scenario['erosion']['TWL']),dhigh)) + 1]
    ax3.set_ylim(ylims)
    
    
    # # Subplot for Morphology Change
    # ax4 = axs[1,2]
    cmap = plt.colormaps['viridis'].resampled(lutsize=len(scenario['erosion']['zmat_times']))
    colors = cmap(np.linspace(0, 1, len(scenario['erosion']['zmat_times'])+1))
    cumDV_erosion = -np.cumsum(scenario['erosion']['dV'])
    cumDV_accretion = np.cumsum(scenario['accretion']['dV'])
    cumDV_net = cumDV_accretion + cumDV_erosion
    
    ylims4 = [np.min(cumDV_erosion)-0.5, np.max(cumDV_accretion)+0.5]
    
    ax4.plot(times, cumDV_erosion, 'b-',linewidth=3, label='Erosion (Waves)')
    ax4.plot(times, cumDV_accretion, 'r-',linewidth=3, label='Accretion (Winds)')
    ax4.plot(times, cumDV_net, 'k--', linewidth=4, label='Net')
    
    for i, t in enumerate(scenario['erosion']['zmat_times']):
        ax4.plot([t, t], ylims4, '--', linewidth=2, color=colors[i])
        
    ax4.set_xlim(xlims)
    ax4.set_ylim(ylims4)
    ax4.set_ylabel(r'$\Delta V_{dune}\ [m^3/m]$')
    ax4.set_title('Dune Volume Change',fontweight='bold')
    ax4.legend(loc='best')
    ax4.grid(True)
    ax4.spines[:].set_linewidth(1.5)
    ax4.xaxis.set_major_formatter(date_fmt)
    ax4.tick_params(axis='x',rotation=30)
    
    # # Subplot for profile change
    # ax5 = axs[1,0:2]
    for i, t in enumerate(scenario['erosion']['zmat_times']):
        ax5.plot(scenario['grids']['XGrid'],scenario['erosion']['Z'][i],'-',color=colors[i],linewidth=3)
    ylims5 = [0, np.nanmax(np.nanmax(scenario['erosion']['Z']))+1]
    ax5.set_ylim(ylims5)
    ax5.set_ylabel('Z [m, NAVD]')
    ax5.set_xlabel('Cross-Shore Disntance [m]')
    ax5.grid(True)
    ax5.set_title('Dune Profile Change (Erosion Only)', fontweight='bold')
    ax5.spines[:].set_linewidth(1.5)
    
    # erosion category:
    if cumDV_net[-1]>-0.5:
        col = [0.2, 0.9, 0.2]
        tex = 'No or Minimal Net Dune Erosion Predicted'
    elif cumDV_net[-1]>-5:
        col = [1, 0, 1]
        tex = 'Minor Dune Erosion Predicted'
    else:
        col = [1, 0, 0]
        tex = 'Substantial Dune erosion predicted'

    # add background patch:
    xgrid = scenario['grids']['XGrid']
    ax5.fill_between(xgrid, ylims5[0], ylims5[1], color=col, alpha=0.075)
    
    if np.min(scenario['grids']['ZGrid'])>0:
        x0 = np.min(scenario['grids']['XGrid'])
    else:
        x0 = np.interp(0,scenario['grids']['XGrid'],scenario['grids']['ZGrid'])
    
    xlim = [x0, np.max(scenario['grids']['XGrid'])]

    # Add additional info onto plot
    xlims5 = xlim
    ax5.set_xlim(xlims5)

    xo = xlims5[0] + (xlims5[1] - xlims5[0])*0.05
    xo2 = xlims5[0] + (xlims5[1] - xlims5[0])*0.5
    yo2 = ylims5[0] + (ylims5[1] - ylims5[0])*0.08
    
    ax5.plot(xgrid, np.full_like(xgrid, dtoe), '--', color=[0.9, 0.3, 0.3], linewidth=2)
    ax5.plot(xgrid, np.full_like(xgrid, dhigh), '--', color=[0.3, 0.1, 0.1], linewidth=2)
    ax5.text(xo, dhigh+0.25, 'Dune Crest', fontweight='bold',color=[0.3, 0.1, 0.1])
    ax5.text(xo, dtoe+0.25, 'Dune Toe', fontweight='bold',color=[0.9, 0.3, 0.3])
    ax5.text(xo2, yo2, tex, fontweight='bold', color=col)
    
    
def drt_frf_plotting(scenario,lidar_dat):
    # drt_plotting: code to plot basic model output for deterministic run  
    fig = plt.figure(figsize=(16, 9))
    gs = gridspec.GridSpec(2, 3, figure=fig, wspace=0.25, hspace=0.35)
    
    # Top row (3 subplots)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    
    # Bottom row:
    ax5 = fig.add_subplot(gs[1, 0:2])  # spans bottom-left two columns
    ax4 = fig.add_subplot(gs[1, 2])    # bottom-right panel
    # figure('units','normalized','outerposition',[0 0 1 1])
    
    # ---- COMMON time limits ----
    times = scenario['timing']['times']
    xlims = [np.min(times), np.max(times)]
    xpts = np.linspace(xlims[0].toordinal(), xlims[1].toordinal(), 5)
    xpts = np.unique(np.floor(xpts))
    
    # Format x-axis for date ticks
    date_fmt = mdates.DateFormatter('%Y-%m-%d')

    # # Subplot for Wind Time Series
    # ax1 = axs[0,0]
    ax1.plot(times, scenario['env']['winds']['windSpeed'],'k-',linewidth=3)
    ax1.set_ylabel('Wind Speed [m/s]')
    ax1.grid(True)
    ax1.set_title('Winds', fontweight='bold')
    ax1.set_xlim(xlims)
    ax1.xaxis.set_major_formatter(date_fmt)
    ax1.tick_params(axis='x',rotation=30)
    ax1.tick_params(width=1.5)
    ax1.spines[:].set_linewidth(1.5)
    
    ax1b = ax1.twinx()
    ax1b.plot(times,scenario['env']['winds']['windDirection'],'.',linewidth=0.5,color=[0.5, 0.5, 0.5],markersize=6) 
    ax1b.set_ylabel(r'Wind Direction $[^{\circ}]$')
    ax1b.set_ylim([0,360])
    ax1.yaxis.label.set_color('k')
    ax1b.yaxis.label.set_color([0.5, 0.5, 0.5])
    
    # # Subplot for Wave Time Series   
    # ax2 = axs[0,1]
    ax2.plot(times, scenario['env']['waves']['Hs_deepwater'],'k-',linewidth=3)
    ax2.set_ylabel(r'$H_s$ [m]')
    ax2.grid(True)
    ax2.set_title('Waves',fontweight='bold')
    ax2.set_xlim(xlims)
    ax2.xaxis.set_major_formatter(date_fmt)
    ax2.tick_params(axis='x',rotation=30)
    ax2.spines[:].set_linewidth(1.5)
    
    ax2b  = ax2.twinx()
    ax2b.plot(times,scenario['env']['waves']['D_deepwater'],'.',color=[0.5, 0.5, 0.5],markersize=6)
    ax2b.set_ylabel(r'Wave Direction $[^{\circ}]$')
    ax2b.set_ylim([0,360])
    ax2.yaxis.label.set_color('k')
    ax2b.yaxis.label.set_color([0.5, 0.5, 0.5])
    
    # # Subplot for Water Level Time Series
    # ax3 = axs[0,2]
    ax3.plot(times, scenario['env']['tides']['wl'],color=[0.5, 0.5, 0.5], linewidth=3,label='SWL')
    ax3.plot(times,scenario['erosion']['TWL'],'k-',linewidth=3,label='TWL')
    
    dtoe = scenario['grids']['morphometrics']['dtoe']
    dhigh = scenario['grids']['morphometrics']['dhigh']
    ax3.plot(xlims,[dtoe, dtoe],'--',color=[0.9, 0.3, 0.3],linewidth=2)
    ax3.plot(xlims,[dhigh, dhigh],'--',color=[0.3, 0.1, 0.1],linewidth=2)
    xo = xlims[0] + (xlims[1]-xlims[0])*0.05
    ax3.text(xo, dhigh+0.25,'Dune Crest',fontweight='bold',color=[0.3, 0.1, 0.1])
    ax3.text(xo, dtoe+0.25,'Dune Toe',fontweight='bold',color=[0.9, 0.3, 0.3])
    
    ax3.set_ylabel('SWL & TWL [m]')
    ax3.set_title('Water Levels',fontweight='bold')
    ax3.grid(True)
    ax3.set_xlim(xlims)
    ax3.spines[:].set_linewidth(1.5)
    ax3.xaxis.set_major_formatter(date_fmt)
    ax3.tick_params(axis='x',rotation=30)
    ylims = [np.min((np.min(scenario['env']['tides']['wl']),dtoe)) - 1, np.max((np.max(scenario['erosion']['TWL']),dhigh)) + 1]
    ax3.set_ylim(ylims)
    
    
    # # Subplot for Morphology Change
    # ax4 = axs[1,2]
    cmap = plt.colormaps['viridis'].resampled(lutsize=len(scenario['erosion']['zmat_times']))
    colors = cmap(np.linspace(0, 1, len(scenario['erosion']['zmat_times'])+1))
    cumDV_erosion = -np.cumsum(scenario['erosion']['dV'])
    cumDV_accretion = np.cumsum(scenario['accretion']['dV'])
    cumDV_net = cumDV_accretion + cumDV_erosion
    
    ylims4 = [np.min(cumDV_erosion)-0.5, np.max(cumDV_accretion)+0.5]
    
    ax4.plot(times, cumDV_erosion, 'b-',linewidth=3, label='Erosion (Waves)')
    ax4.plot(times, cumDV_accretion, 'r-',linewidth=3, label='Accretion (Winds)')
    ax4.plot(times, cumDV_net, 'k--', linewidth=4, label='Net')
    
    for i, t in enumerate(scenario['erosion']['zmat_times']):
        ax4.plot([t, t], ylims4, '--', linewidth=2, color=colors[i])
        
    ax4.set_xlim(xlims)
    ax4.set_ylim(ylims4)
    ax4.set_ylabel(r'$\Delta V_{dune}\ [m^3/m]$')
    ax4.set_title('Dune Volume Change',fontweight='bold')
    ax4.legend(loc='best')
    ax4.grid(True)
    ax4.spines[:].set_linewidth(1.5)
    ax4.xaxis.set_major_formatter(date_fmt)
    ax4.tick_params(axis='x',rotation=30)
    
    # # Subplot for profile change
    # ax5 = axs[1,0:2]
    for i, t in enumerate(scenario['erosion']['zmat_times']):
        ax5.plot(scenario['grids']['XGrid'],scenario['erosion']['Z'][i],'-',color=colors[i],linewidth=3)
    ylims5 = [0, np.nanmax(np.nanmax(scenario['erosion']['Z']))+1]
    ax5.set_ylim(ylims5)
    ax5.set_ylabel('Z [m, NAVD]')
    ax5.set_xlabel('Cross-Shore Disntance [m]')
    ax5.grid(True)
    ax5.set_title('Dune Profile Change (Erosion Only)', fontweight='bold')
    ax5.spines[:].set_linewidth(1.5)
    
    # erosion category:
    if cumDV_net[-1]>-0.5:
        col = [0.2, 0.9, 0.2]
        tex = 'No or Minimal Net Dune Erosion Predicted'
    elif cumDV_net[-1]>-5:
        col = [1, 0, 1]
        tex = 'Minor Dune Erosion Predicted'
    else:
        col = [1, 0, 0]
        tex = 'Substantial Dune erosion predicted'

    # add background patch:
    xgrid = scenario['grids']['XGrid']
    ax5.fill_between(xgrid, ylims5[0], ylims5[1], color=col, alpha=0.075)
    
    if np.min(scenario['grids']['ZGrid'])>0:
        x0 = np.min(scenario['grids']['XGrid'])
    else:
        x0 = np.interp(0,scenario['grids']['XGrid'],scenario['grids']['ZGrid'])
    
    xlim = [x0, np.max(scenario['grids']['XGrid'])]

    # Add additional info onto plot
    xlims5 = xlim
    ax5.set_xlim(xlims5)

    xo = xlims5[0] + (xlims5[1] - xlims5[0])*0.05
    xo2 = xlims5[0] + (xlims5[1] - xlims5[0])*0.5
    yo2 = ylims5[0] + (ylims5[1] - ylims5[0])*0.08
    
    ax5.plot(xgrid, np.full_like(xgrid, dtoe), '--', color=[0.9, 0.3, 0.3], linewidth=2)
    ax5.plot(xgrid, np.full_like(xgrid, dhigh), '--', color=[0.3, 0.1, 0.1], linewidth=2)
    ax5.text(xo, dhigh+0.25, 'Dune Crest', fontweight='bold',color=[0.3, 0.1, 0.1])
    ax5.text(xo, dtoe+0.25, 'Dune Toe', fontweight='bold',color=[0.9, 0.3, 0.3])
    ax5.text(xo2, yo2, tex, fontweight='bold', color=col)
