# coding: utf-8
from __future__ import print_function

import math
import numpy as np
import os
import pandas as pd
import netCDF4 as nc
import logging
from scipy import interpolate

import dassflow1d.m_obs as m_obs


class RiverNodesOverpass:
    
    def __init__(self, X, H, W):
        self.X = X
        self.H = H
        self.W = W


def sge_rivernodes_observations(passplanfile, datadir, prefix, mesh):
    """ Create observation from SGE RiverNodes files
    
        :param X: Curvilinear abscissae of the cross-sections
        :type X: numpy.ndarray ([x])
        :param B: Bathymetry elevations of the cross-sections
        :type B: numpy.ndarray ([x])
        :param Z: Observed heights of the cross-sections
        :type Z: numpy.ndarray ([t, x])
        :param W: Observed widths of the cross-sections
        :type W: numpy.ndarray ([t, x])
        :param dx: dx Spacing between cross-sections
        :type dx: float
        :param dZmin: Threshold for cleaning close observations
        :type dZmin: float
        
        :return Observations object
        
    """
  
    # CHECK-UP
    # TODO
    
    # Load passplan
    passplan = pd.read_csv(passplanfile, delim_whitespace=True, header=4)

    # Compute simulation times
    t = (passplan["DayOfYear"] - 1.0) * 86400.0
    
    # Load first data file to get number of nodes
    cycle = passplan["cycle"][0]
    orbit = passplan["orbit"][0]
    day = int(passplan["DayOfYear"][0])
    data = __read_nodes_file__(datadir, prefix, cycle, orbit, day)
    
    # Initialise Observations object
    obs = m_obs.Observations(data.X.size)
    
    # Setup stations
    __setup_stations__(obs, mesh, t, data.X)
      
    # Setup observations data
    __setup_observations_data__(obs, datadir, prefix, passplan)
    
    return obs

    
def __read_nodes_file__(datadir, prefix, cycle, orbit, day):
    """ Read a RiverNodes data file
    
        TODO : document arguments
    """
        
    fname = os.path.join(datadir, "%s-Pass-%i-Day-%i.nc" % (prefix, orbit, int(day)))
    
    dataset = nc.Dataset(fname)
    
    # Retrieve curvilinear abscissae of nodes
    X = dataset.variables["xs"][:]
    H = dataset.variables["height"][:]
    W = dataset.variables["width"][:]
    
    # Create RiverNodes overpass
    data = RiverNodesOverpass(X, H, W)
    
    return data

    
def __setup_stations__(obs, mesh, t, X, dxmax=None):
    """ Setup cross-sections
    
        TODO : document arguments
    """
    
    xcs = np.zeros(mesh.ncs)
    for ics in range(0, mesh.ncs):
      xcs[ics] = mesh.cs[ics].x

    sta_cs = []
    for ista in range(0, X.size):
      sta_cs.append([])

    Xs = X[-1] - X
      
    # Associate mesh cross-sections to stations
    ista = np.ones(mesh.ncs, dtype=int) * -1
    for ics in range(0, mesh.ncs):
      dx = np.abs(mesh.cs[ics].x - Xs)
      imin = np.argmin(dx)
      if dxmax is not None:
        if dx[imin] < dxmax:
          ista[ics] = imin
          sta_cs[imin].append(ics)
      else:
        ista[ics] = imin
        sta_cs[imin].append(ics)
        
    # Setup stations
    for ista in range(0, X.size):
      obs.stations[ista].setup_station(sta_cs[ista], t)
      obs.stations[ista].offset = ista * t.size

    
def __setup_observations_data__(obs, datadir, prefix, passplan):
    """ Setup observations data
    
        TODO : document arguments
    """
    
    # Compute number of data
    ndata = len(obs.stations) * passplan.shape[0]
    obs.setup_observations_data(ndata)

    
    # Loop on passes
    for ipass in range(0, passplan.shape[0]):
        
        # Retrieve cycle, orbit and day
        cycle = passplan["cycle"][ipass]
        orbit = passplan["orbit"][ipass]
        day = int(passplan["DayOfYear"][ipass])
        
        # Load data file
        data = __read_nodes_file__(datadir, prefix, cycle, orbit, day)
        
        # Check same number of stations in every overpasses
        if data.X.size != len(obs.stations):
            raise RuntimeError("Overpass with different number of stations")
        
        # Update observations
        print(obs.obs.shape)
        for ista in range(0, data.X.size):
            obs.obs[0, obs.stations[ista].offset+ipass] = data.H[ista]
            obs.obs[1, obs.stations[ista].offset+ipass] = data.W[ista]
