# coding: utf-8
from __future__ import print_function

import logging
import netCDF4 as nc
import numpy as np
import os


from .SwotObs import SwotObs


class PepsiNetCDF(SwotObs):
  """ Class for parsing PEPSI data in NetCDF format
  """
  
  def __init__(self, fname, obs_only=True, reach_data=False):
    """ Instanciate a PepsiNetCDF object
    
        :param fname: Path to the data file
        :type fname: str
        :param obs_only: If True only observations are read
        :type obs_only: bool
        :param reach_data: If True observations at reach scale are read
        :type reach_data: bool
        
    """
    
    # Retrieve logger and append debug messages
    logger = logging.getLogger("HiVDI")
    logger.debug("Instanciate PepsiNetCDF object <%s>" % id(self))
    logger.debug("- fname: %s" % fname)
  
    # CHECK-UP
    if not isinstance(fname, str):
      raise ValueError("'fname' must be a string")
    if not os.path.isfile(fname):
      raise IOError("File not found: %s" % fname)

    # Call to ancestor's __init__
    super(PepsiNetCDF, self).__init__()
      
    # Open NetCDF file
    dataset = nc.Dataset(fname)

    # Retrieve 'River_Info' group
    group = dataset.groups['River_Info']

    # Retrieve QWBM
    self.QWBM = group.variables['QWBM'][0]
    
    # Read reach data
    if reach_data:
      
      # Compute curvilinear abscissae
      rch_bnd = group.variables['rch_bnd'][:]
      self.good = group.variables['gdrch'][:].astype(dtype=int)
      self.X = 0.5 * (rch_bnd[0:-1] + rch_bnd[1:])

      # Retrieve 'Reach_Timeseries' group
      group = dataset.groups['Reach_Timeseries']
      
      # Parse observation variables
      self.t = self.__parse_t_variable__("t", group)
      self.t0 = self.t[0]
      self.Z = self.__parse_xt_variable__("H", group)
      self.W = self.__parse_xt_variable__("W", group)
      self.S = self.__parse_xt_variable__("S", group)
      
      # Parse extra variables
      if not obs_only:
        self.A = self.__parse_xt_variable__("A", group)
        self.P = self.__parse_xt_variable__("P", group)
        self.Q = self.__parse_xt_variable__("Q", group)
    
    # Read cross-sections data
    else:
      
      # Retrieve 'XS_Timeseries' group
      group = dataset.groups['XS_Timeseries']
      self.xs_rch = self.__parse_x_variable__("xs_rch", group)
      
      # Parse observation variables
      self.t = self.__parse_t_variable__("t", group)
      self.t0 = self.t[0]
      self.X = self.__parse_x_variable__("X", group)
      self.Z = self.__parse_xt_variable__("H", group)
      self.W = self.__parse_xt_variable__("W", group)
      
      ## Compute S
      reaches_group = dataset.groups['Reach_Timeseries']
      Sr = self.__parse_xt_variable__("S", reaches_group)
      thresholds = np.zeros(self.Z.shape)
      for ix in range(0, self.X.size):
        ir = self.xs_rch[ix]
        thresholds[:, ix] = 0.5 * Sr[:, ir-1]
      SwotObs.compute_S(self, None, thresholds, None)
      
      # Parse extra variables
      if not obs_only:
        self.B = self.__parse_x_variable__("Z", group)
        self.Q = self.__parse_xt_variable__("Q", group)
        self.A = self.__parse_xt_variable__("A", group)
        self.P = self.__parse_xt_variable__("P", group)
        self.n = self.__parse_xt_variable__("n", group)
        
    self.H = self.Z

    
  def compute_S(self, filtering=None, threshold=1e-9, validity_flag=False):
    
    if threshold == "auto":
      if self.reaches is None:
        raise RuntimeError("Cannot compute automatic slope thresholds without reaches data")
      
      thresholds = np.zeros((self.nt, self.nx))
      for ix in range(0, self.nx):
        ir = self.xs_rch[ix]
        thresholds[:, ix] = 0.5 * self.reaches.S[:, ir-1]
      threshold = thresholds
    
    ModelResults.compute_S(self, filtering, threshold, validity_flag)
    #print("self.S=", self.S)
    
    if self.reaches is not None:
      
      if self.reaches.S is not None and isinstance(threshold, float):
        self.reaches.S = np.maximum(self.reaches.S, threshold)


  def __parse_t_variable__(self, varname, dataset):
    """ Parse a time varying variable. Handles masked array with replacement 
        of masked values with numpy.nan
    
        :param varname: Name of the variable in dataset
        :type varname: str
        :param dataset: dataset object
        :type dataset: netCDF4.DataSet
    """
    
    array = dataset.variables[varname][:, 0] * 86400.0
    if isinstance(array, np.ma.core.MaskedArray):
      array = array.filled(fill_value=np.nan)
    
    return array


  def __parse_x_variable__(self, varname, dataset):
    """ Parse a space varying variable. Handles masked array with replacement 
        of masked values with numpy.nan
    
        :param varname: Name of the variable in dataset
        :type varname: str
        :param dataset: dataset object
        :type dataset: netCDF4.DataSet
    """
    
    #array = dataset.variables[varname][0, ::-1]
    array = dataset.variables[varname][0, :]
    if isinstance(array, np.ma.core.MaskedArray):
      array = array.filled(fill_value=np.nan)
    
    return array


  def __parse_xt_variable__(self, varname, dataset):
    """ Parse a space and time varying variable. Handles masked array with 
        replacement of masked values with numpy.nan
    
        :param varname: Name of the variable in dataset
        :type varname: str
        :param dataset: dataset object
        :type dataset: netCDF4.DataSet
    """
    
    if not varname in dataset.variables:
      return None
    
    array = dataset.variables[varname][:, :]
    if isinstance(array, np.ma.core.MaskedArray):
      array = array.filled(fill_value=np.nan)
    
    return array
