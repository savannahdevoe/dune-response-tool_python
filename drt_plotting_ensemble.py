# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 12:03:38 2025

@author: sdevo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.cm import get_cmap
from datetime import datetime

def drt_plotting_ensemble(scenario, scenarioLow, scenarioHigh):
    # drt_plotting_enseble(scenario,senarioLow,scenarioHigh):

    # Convert to numpy datetime64 if needed
    times = np.array(scenario['timing']['times'], dtype='datetime64[ns]')
    xlims = [times.min(), times.max()]
    xpts = np.linspace(xlims[0], xlims[1], 5)

    fig, axs = plt.subplots(2, 3, figsize=(14, 8))
    fig.subplots_adjust(wspace=0.3, hspace=0.4)
    
    # === (1) Winds ===
    ax = axs[0, 0]
    ax.plot(times, scenario['env']['winds']['windSpeed'], 'k-', linewidth=3)
    ax.set_xlim(xlims)
    ax.set_ylabel('Wind Speed (m/s)')
    ax.grid(True, linewidth=1.5)
    ax.set_title('Winds', fontweight='bold')

    ax2 = ax.twinx()
    ax2.plot(times, scenario['env']['winds']['windDirection'], '.', color=[0.5, 0.5, 0.5], markersize=8)
    ax2.set_ylabel('Wind Direction (°)')
    ax2.set_ylim([0, 360])

    ax.tick_params(axis='x', rotation=30)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))

    # === (2) Waves ===
    ax = axs[0, 1]
    ax.plot(times, scenario['env']['waves']['Hs_deepwater'], 'k-', linewidth=3)
    ax.set_xlim(xlims)
    ax.set_ylabel('Hₛ (m)')
    ax.grid(True, linewidth=1.5)
    ax.set_title('Waves', fontweight='bold')

    ax2 = ax.twinx()
    ax2.plot(times, scenario['env']['waves']['D_deepwater'], '.', color=[0.5, 0.5, 0.5], markersize=8)
    ax2.set_ylabel('Wave Direction (°)')
    ax2.set_ylim([0, 360])
    ax.tick_params(axis='x', rotation=30)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d'))

    # === (3) Water Levels ===
    ax = axs[0, 2]
    ax.plot(times, scenario['env']['tides']['wl'], color=[0.5, 0.5, 0.8], linewidth=3, label='SWL')
    ax.plot(times, scenario['erosion']['TWL'], 'k-', linewidth=3, label='TWL')
    ax.plot(times, scenarioLow['erosion']['TWL'], 'k-', linewidth=1)
    ax.plot(times, scenarioHigh['erosion']['TWL'], 'k-', linewidth=1)

    # Ensemble shading
    max_twl = np.nanmax(np.vstack([
        scenario['erosion']['TWL'],
        scenarioLow['erosion']['TWL'],
        scenarioHigh['erosion']['TWL']
    ]), axis=0)
    min_twl = np.nanmin(np.vstack([
        scenario['erosion']['TWL'],
        scenarioLow['erosion']['TWL'],
        scenarioHigh['erosion']['TWL']
    ]), axis=0)

    ax.fill_between(times, min_twl, max_twl, color='k', alpha=0.1, label='Ensemble TWL')

    dtoe = scenario['grids']['morphometrics']['dtoe']
    dhigh = scenario['grids']['morphometrics']['dhigh']
    ax.axhline(dtoe, ls='--', color=[0.9, 0.3, 0.3], linewidth=2)
    ax.axhline(dhigh, ls='--', color=[0.3, 0.1, 0.1], linewidth=2)
    ax.text(times[0], dhigh + 0.25, 'Dune Crest', fontweight='bold', color=[0.3, 0.1, 0.1])
    ax.text(times[0], dtoe + 0.25, 'Dune Toe', fontweight='bold', color=[0.9, 0.3, 0.3])

    ax.set_ylabel('SWL & TWL (m)')
    ax.set_xlim(xlims)
    ax.grid(True)
    ax.set_title('Water Levels', fontweight='bold')
    ax.legend(loc='upper right')

    # === (4) Dune Volume Change ===
    ax = axs[1, 2]
    cmap = get_cmap('viridis', len(scenario['erosion']['zmat_times']) + 1)

    def cum_volumes(scn):
        erosion = -np.cumsum(scn['erosion']['dV'])
        accretion = np.cumsum(scn['accretion']['dV'])
        return erosion, accretion, erosion + accretion

    ero, acc, net = cum_volumes(scenario)
    ero_low, acc_low, net_low = cum_volumes(scenarioLow)
    ero_high, acc_high, net_high = cum_volumes(scenarioHigh)

    max_net = np.nanmax(np.vstack([net, net_low, net_high]), axis=0)
    min_net = np.nanmin(np.vstack([net, net_low, net_high]), axis=0)
    ax.fill_between(times, min_net, max_net, color='k', alpha=0.1, label='Ensemble Envelope')
    ax.plot(times, net, 'k-', linewidth=3)

    for itime, zt in enumerate(scenario['erosion']['zmat_times']):
        ax.axvline(zt, ls='--', linewidth=2, color=cmap(itime))

    ax.set_ylabel('ΔV₍dune₎ (m³/m)')
    ax.set_xlim(xlims)
    ax.grid(True)
    ax.set_title('Dune Volume Change', fontweight='bold')
    ax.legend(loc='lower left')

    # === (5) Dune Profile Change ===
    ax = axs[1, 0]
    for itime, zt in enumerate(scenario['erosion']['zmat_times']):
        ax.plot(scenario['grids']['XGrid'], scenario['erosion']['Z'][:, itime], '-', color=cmap(itime), linewidth=3)
    ax.set_ylabel('Z (m, NAVD)')
    ax.set_xlabel('Cross-Shore Distance (m)')
    ax.grid(True)
    ax.set_title('Dune Profile Change (Base Case)', fontweight='bold')

    # Determine erosion classification
    if net[-1] > -0.5:
        col, tex = [0.2, 0.9, 0.2], 'No or Minimal Net Dune Erosion Predicted'
    elif net[-1] > -5:
        col, tex = [1, 0, 1], 'Minor Dune Erosion Predicted'
    else:
        col, tex = [1, 0, 0], 'Substantial Dune Erosion Predicted'

    ylims = ax.get_ylim()
    xlims = ax.get_xlim()
    ax.fill_betweenx(ylims, xlims[0], xlims[1], color=col, alpha=0.075)
    ax.text(xlims[0] + 0.5*(xlims[1]-xlims[0]), ylims[0] + 0.1*(ylims[1]-ylims[0]), tex, fontweight='bold', color=col)

    plt.show()