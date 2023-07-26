# coding: utf-8
from __future__ import print_function

import math
import numpy as np
import logging
from scipy import interpolate

import dassflow1d.m_obs as m_obs


def pepsi_nodes_observations(pepsidata, mesh):
    """ Create observation from PEPSI data
    
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
    
    # Initialise Observations object
    obs = m_obs.Observations(pepsidata.X.size)
    
    # Setup stations
    __setup_stations__(obs, mesh, pepsidata.t, pepsidata.X)
      
    # Setup observations data
    __setup_observations_data__(obs, pepsidata.t, pepsidata.X, pepsidata.H, pepsidata.W)
    
    return obs

    
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

    
def __setup_observations_data__(obs, t, X, H, W):
    """ Setup observations data
    
        TODO : document arguments
    """
    
    # Compute number of data
    ndata = X.size * t.size
    
    # Allocate observations data
    obs.setup_observations_data(ndata)
    
    # Compute data
    idata = 0
    for ista in range(0, X.size):
      Hx = H[:, ista]
      Wx = W[:, ista]
      obs.stations[ista].offset = idata
      obs.obs[0, idata:idata+t.size] = Hx
      obs.obs[1, idata:idata+t.size] = Wx
      idata += t.size
