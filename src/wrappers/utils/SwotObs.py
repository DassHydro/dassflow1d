# coding: utf-8
from __future__ import print_function

import numpy as np
import logging


class SwotObs(object):
  """ Class for handling SWOT observations
  """
  
  def __init__(self, t=None, X=None, W=None, Z=None):
    """ Instanciate a SwotObs object
    
        :param t: Datetimes
        :type t: array or list
        :param X: Curvilinear abscissae
        :type X: array
        :param W: Observed widths (rows:time, cols:space)
        :type W: array
        :param Z: Observed elevations (rows:time, cols:space)
        :type Z: array
        
    """
    
    # Retrieve logger and append debug messages
    logger = logging.getLogger("HiVDI")
    logger.debug("Instanciate SwotObs object <%s>" % id(self))
    logger.debug("- X=%s" % str(X))
    logger.debug("- W=%s" % str(W))
    logger.debug("- Z=%s" % str(Z))
  
    # CHECK-UP
    if t is not None:
      if isinstance(t, list):
        dates = t
        if isinstance(t[0], str):
          t = np.zeros(len(dates), dtype='datetime64[s]')
          for it in range(0, len(dates)):
            t[it] = np.datetime64(dates[it]).astype('datetime64[s]')
        else:
          t = np.array(dates)
      elif not isinstance(t, np.ndarray):
        try:
          t = np.array(t)
          t = t.astype("float")
          if t.ndim == 0:
            t = t.reshape(1)
        except:
          raise ValueError("'t' must be a 1d array")
      if isinstance(t, np.ndarray):
        if t.ndim != 1:
          raise ValueError("'t' must be a 1d array")
    if X is not None:
      if not isinstance(X, np.ndarray):
        try:
          X = np.array(X)
          X = X.astype("float")
          if X.ndim == 0:
            X = X.reshape(1)
        except:
          raise ValueError("'X' must be a 1d array")
      if isinstance(X, np.ndarray):
        if X.ndim != 1:
          raise ValueError("'X' must be a 1d array")
    if X is not None and W is None:
      raise ValueError("'W' must be specified if 'X' is specified")
    if W is not None:
      if not isinstance(W, np.ndarray):
        try:
          W = np.array(W)
          W = W.astype("float")
          if W.ndim == 0:
            W = W.reshape((1,1))
        except:
          raise ValueError("'W' must be a 2d array")
      if isinstance(W, np.ndarray):
        if W.ndim != 2:
          raise ValueError("'W' must be a 2d array")
        if W.shape[0] != t.size or W.shape[1] != X.size:
          raise ValueError("'W' must be a of shape (number of dates, number "
                           "of points)")
    if X is not None and Z is None:
      raise ValueError("'Z' must be specified if 'X' is specified")
    if Z is not None:
      if not isinstance(Z, np.ndarray):
        try:
          Z = np.array(Z)
          Z = Z.astype("float")
          if Z.ndim == 0:
            Z = Z.reshape((1,1))
        except:
          raise ValueError("'Z' must be a 2d array")
      if isinstance(Z, np.ndarray):
        if Z.ndim != 2:
          raise ValueError("'Z' must be a 2d array")
        if Z.shape[0] != t.size or Z.shape[1] != X.size:
          raise ValueError("'Z' must be a of shape (number of dates, number "
                           "of points)")
      
    # Set all properties
    self.t = t
    self.X = X
    self.W = W
    self.Z = Z
    
    
  def compute_dA(self):
    """ Compute changes in flow area (dA)
    """
    
    # Retrieve logger and append debug messages
    logger = logging.getLogger("HiVDI")
    logger.debug("Compute dA for SwotObs object <%s>" % id(self))
    
    # Special treatment for empty observations set
    if self.W is None or self.Z is None:
      self.dA = None
      return
    
    # Init array of changes in flow area
    self.dA = np.zeros(self.Z.shape)
    
    # Loop on points
    for ix in range(0, self.Z.shape[1]):
      
      # Sort times by increasing elevations
      sorted_t = np.argsort(self.Z[:, ix])
      
      for it in range(1, sorted_t.size):
      
        Wmean = 0.5 * (self.W[sorted_t[it], ix] + self.W[sorted_t[it-1], ix])
        dZ = self.Z[sorted_t[it], ix] - self.Z[sorted_t[it-1], ix]
        self.dA[sorted_t[it], ix] = self.dA[sorted_t[it-1], ix] + Wmean * dZ
    
    
  def compute_A(self):
    """ Compute flow areas
    """
    
    # Retrieve logger and append debug messages
    logger = logging.getLogger("HiVDI")
    logger.debug("Compute A for SwotObs object <%s>" % id(self))
    
    # CHECK-UP
    if not hasattr(self, "B"):
      if not hasattr(self, "A0"):
        raise ValueError("SwotObs object must have A0 or B to compute flow areas")
    if not hasattr(self, "dA"):
      raise ValueError("SwotObs object must have dA to compute flow areas")
    
    # Compute unobserved flow area if only B attribute is set
    if hasattr(self, "B"):
      self.compute_A0()
      
    # Compute flow areas
    self.A = self.dA.copy()
    for ix in range(0, self.X.size):
      self.A[:, ix] += self.A0[ix]
    
    
  def compute_A0(self):
    """ Compute unobserved flow areas
    """
    
    # Retrieve logger and append debug messages
    logger = logging.getLogger("HiVDI")
    logger.debug("Compute A0 for SwotObs object <%s>" % id(self))
    
    # CHECK-UP
    if not hasattr(self, "B"):
      if not hasattr(self, "A"):
        raise ValueError("SwotObs object must have A or B to compute A0")
      elif not hasattr(self, "dA"):
        raise ValueError("SwotObs object must have dA to compute A0 from A")
    
    # Compute unobserved flow area if B attribute is set
    if hasattr(self, "B"):
      Z0 = np.nanmin(self.Z, axis=0)
      W0 = np.nanmin(self.W, axis=0)
      self.A0 = (Z0 - self.B) * W0
      
    # Compute unobserved flow area using A and dA
    else:
      self.A0 = self.A[0, :] - self.dA[0, :]
    
    
  def compute_B(self):
    """ Compute bathymetry
    """
    
    # Retrieve logger and append debug messages
    logger = logging.getLogger("HiVDI")
    logger.debug("Compute B for SwotObs object <%s>" % id(self))
    
    # CHECK-UP
    if not hasattr(self, "A0"):
      if not hasattr(self, "A"):
        raise ValueError("SwotObs object must have A or A0 to compute B")
    
    # Compute bathymetry if A0 attribute is set
    if hasattr(self, "A0"):
      Z0 = np.nanmin(self.Z, axis=0)
      W0 = np.nanmin(self.W, axis=0)
      self.B = Z0 - self.A0 / W0
      
    # Otherwise compute bathymetry using A
    else:
      Z0 = np.nanmin(self.Z, axis=0)
      W0 = np.nanmin(self.W, axis=0)
      A0 = np.nanmin(self.A, axis=0)
      self.B = Z0 - A0 / W0

      
  def compute_S(self, filtering=None, threshold=1e-9, validity_flag=False):
    """ Compute slopes
    """
    
    # Retrieve logger and append debug messages
    logger = logging.getLogger("HiVDI")
    logger.debug("Compute S for SwotObs object <%s>" % id(self))
    
    self.S = np.ones(self.Z.shape, dtype=float)
    if validity_flag:
      self.Svalid = np.ones(self.Z.shape, dtype=bool)

    if isinstance(threshold, np.ndarray):
      thresholds = threshold
    else:
      thresholds = None
      
    #--------------------------------------------------------------------------
    # Loop on x
    #--------------------------------------------------------------------------
    for ix in range(0, self.X.size):
      
      if ix == 0:
        x1 = self.X[ix]
        z1 = self.Z[:, ix]
        x2 = self.X[ix+1]
        z2 = self.Z[:, ix+1]
      #elif ix == self.nx - 1:
        #x1 = self.X[ix-1]
        #z1 = self.Z[:, ix-1]
        #x2 = self.X[ix]
        #z2 = self.Z[:, ix]
      else:
        #x1 = self.X[ix-1]
        #z1 = self.Z[:, ix-1]
        #x2 = self.X[ix+1]
        #z2 = self.Z[:, ix+1]
        x1 = self.X[ix-1]
        z1 = self.Z[:, ix-1]
        x2 = self.X[ix]
        z2 = self.Z[:, ix]
      if thresholds is not None:
        self.S[:, ix] = np.maximum(thresholds[:, ix], (z1 - z2) / abs(x2 - x1))
      else:
        self.S[:, ix] = np.maximum(threshold, (z1 - z2) / abs(x2 - x1))
        
      if validity_flag:
        self.Svalid[:, ix] = ((z1 - z2) / abs(x2 - x1) >= threshold)
    
    
  def select_spatial_subset(self, subset):
    """ Compute bathymetry
    """
    
    self.X = self.X[subset[0]:subset[1]+1]
    self.W = self.W[:, subset[0]:subset[1]+1]
    self.Z = self.Z[:, subset[0]:subset[1]+1]
    if hasattr(self, 'A'):
      if self.A is not None:
        self.A = self.A[:, subset[0]:subset[1]+1]
    if hasattr(self, 'A0'):
      if self.A0 is not None:
        self.A0 = self.A0[subset[0]:subset[1]+1]
    if hasattr(self, 'B'):
      if self.B is not None:
        self.B = self.B[subset[0]:subset[1]+1]
    if hasattr(self, 'dA'):
      if self.dA is not None:
        self.dA = self.dA[:, subset[0]:subset[1]+1]
    if hasattr(self, 'Q'):
      if self.Q is not None:
        self.Q = self.Q[:, subset[0]:subset[1]+1]
    if hasattr(self, 'S'):
      if self.S is not None:
        self.S = self.S[:, subset[0]:subset[1]+1]
    
    
  def select_time_window(self, window):
    """ Compute bathymetry
    """
    
    self.t = self.t[window[0]:window[1]+1]
    self.W = self.W[window[0]:window[1]+1, :]
    self.Z = self.Z[window[0]:window[1]+1, :]
    if hasattr(self, 'A'):
      if self.A is not None:
        self.A = self.A[window[0]:window[1]+1, :]
    if hasattr(self, 'dA'):
      if self.dA is not None:
        self.dA = self.dA[window[0]:window[1]+1, :]
    if hasattr(self, 'Q'):
      if self.Q is not None:
        self.Q = self.Q[window[0]:window[1]+1, :]
    if hasattr(self, 'S'):
      if self.S is not None:
        self.S = self.S[window[0]:window[1]+1, :]
    
    
  def select_time_subset(self, subset):
    """ TODO
    """
    
    self.t = self.t[subset[0]:subset[1]+1]
    self.W = self.W[subset[0]:subset[1]+1, :]
    self.Z = self.Z[subset[0]:subset[1]+1, :]
    if hasattr(self, 'A'):
      if self.A is not None:
        self.A = self.A[subset[0]:subset[1]+1, :]
    if hasattr(self, 'dA'):
      if self.dA is not None:
        self.dA = self.dA[subset[0]:subset[1]+1, :]
    if hasattr(self, 'Q'):
      if self.Q is not None:
        self.Q = self.Q[subset[0]:subset[1]+1, :]
    if hasattr(self, 'S'):
      if self.S is not None:
        self.S = self.S[subset[0]:subset[1]+1, :]
    
    
  def dt(self, index=None, t0=None):
    """ Compute timedelta from initial time (internal or argument 't0')
    
        :param index: index of the time occurence
        :type index: int or tuple (window)
        :param t0: custom initial time
        :param t0: float or numpy.datetime64
    """
    
    # Retrieve logger and append debug messages
    #logger = logging.getLogger("HiVDI")
    #logger.debug("Call to SwotObs<%s>::dt" % id(self))
    #logger.debug("-- index=%s" % str(index))
    #logger.debug("-- t0=%s" % str(t0))
    
    # CHECK-UP
    # TODO (check that t0 is set (either internal or custom), etc.
    
    # Compute bathymetry if A0 attribute is set
    if isinstance(index, tuple):
      
      indices = slice(index[0], index[1])
      
    elif index is not None:
      
      indices = index
      
    else:
      
      indices = slice(None)
    
    if self.t.dtype == "datetime64[s]":
      
      # times are in datetime64[s] format
      if t0 is not None:
        
        # CHECK-UP : check that custom t0 is of type numpy.datetime64[s]
        if not isinstance(t0, numpy.datetime64):
          raise ValueError("'t0' must be of type <numpy.datetime64>")
        
        # Return timedelta in seconds since custom t0
        return (self.t[indices] - t0) / np.timedelta64(1, "s")
      
      else:
        
        # Return timedelta in seconds since internal t0
        return (self.t[indices] - self.t0) / np.timedelta64(1, "s")
        
    else:
        
      # times are in float format
      if t0 is not None:
        
        # CHECK-UP : check that custom t0 is of type float
        if not isinstance(t0, float):
          raise ValueError("'t0' must be of type <float>")
        
        # Return timedelta in seconds since custom t0
        return self.t[indices] - t0
      
      else:
        
        # Return timedelta in seconds since internal t0
        return self.t[indices] - self.t0

    #else:
      
      #if hasattr(self.t, "dtype"):
        
        #if self.t.dtype == "datetime64[s]":
          
          ## times are in datetime64[s] format
          #if t0 is not None:
            
            ## CHECK-UP : check that custom t0 is of type numpy.datetime64[s]
            #if not isinstance(t0, numpy.datetime64):
              #raise ValueError("'t0' must be of type <numpy.datetime64>")
            
            ## Return timedelta in seconds since custom t0
            #dt = self.t[:] - t0) / np.timedelta64(1, "s")
          
          #else:
            
            ## Return timedelta in seconds since internal t0
            #dt = (self.t[:] - self.t0) / np.timedelta64(1, "s")
          
      #else:
          
        ## times are in datetime64[s] format
        #if t0 is not None:
          
          ## CHECK-UP : check that custom t0 is of type float
          #if not isinstance(t0, float):
            #raise ValueError("'t0' must be of type <float>")
          
          ## Return timedelta in seconds since custom t0
          #return self.t[index] - t0
        
        #else:
          
          ## Return timedelta in seconds since internal t0
          #return self.t[index] - self.t0
      
