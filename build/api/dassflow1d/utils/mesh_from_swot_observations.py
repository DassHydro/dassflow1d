# coding: utf-8
from __future__ import print_function

import math
import numpy as np
import logging
from scipy import interpolate

import dassflow1d.m_mesh as m_mesh


def mesh_from_swot_observations(X, B, Z, W, dx=None, dZmin=0.0, with_ghost=True):
    """ Create a mesh from SWOT observations
    
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
        
    """
  
    # CHECK-UP
    # TODO
    
    # Initialise mesh
    if dx is None:
        if with_ghost:
            mesh = m_mesh.Mesh(X.size+4, 1)
        else:
            mesh = m_mesh.Mesh(X.size, 1)
    else:
        L = X[-1] - X[0]
        nreg = int(math.ceil(L/dx))
        Lreg = dx * nreg
        Xreg = np.arange(X[0], X[0]+Lreg+dx, dx)
        if with_ghost:
            mesh = m_mesh.Mesh(Xreg.size+4, 1)
        else:
            mesh = m_mesh.Mesh(Xreg.size, 1)
    
    # Create cross-sections
    if dx is None:
        __setup_cross_sections__(mesh, X, B, Z, W, dZmin, with_ghost)
    else:
        __setup_resampled_cross_sections__(mesh, X, B, Z, W, dx, dZmin, with_ghost)
      
    __finalize_mesh__(mesh, with_ghost)
    
    if with_ghost:
        mesh.setup_segment(0, first_cs=2, ncs=mesh.ncs-4, us_seg=[-1], ds_seg=-1) 
    else:
        mesh.setup_segment(0, first_cs=0, ncs=mesh.ncs, us_seg=[-1], ds_seg=-1) 
    
    
    return mesh

    
def __setup_cross_sections__(mesh, X, B, Z, W, dZmin, with_ghost):
    """ Setup cross-sections
    
        TODO : document arguments
    """
    
    if with_ghost:
        offset = 2
    else:
        offset = 0
    
    for ix in range(0, X.size):

        # Retrieve observed Z and W for current cross-section
        Zx = Z[:, ix]
        Wx = W[:, ix]
        
        # Sort observed Z and W
        isort = np.argsort(Zx)
        Zx = Zx[isort]
        Wx = Wx[isort]
        
        # Remove doublons using dZmin threshold
        nodoublons = Zx[1:] > Zx[0:-1] + dZmin
        tmp = np.zeros(np.sum(nodoublons)+1)
        tmp[0] = Zx[0]
        tmp[1:] = (Zx[1:])[nodoublons]
        Zx = tmp
        tmp = np.zeros(np.sum(nodoublons)+1)
        tmp[0] = Wx[0]
        tmp[1:] = (Wx[1:])[nodoublons]
        Wx = tmp
        
        mesh.cs[ix+offset].x = X[-1] - X[ix]
        mesh.cs[ix+offset].coord.x = X[-1] - X[ix]
        mesh.cs[ix+offset].coord.y = 0.0
        mesh.cs[ix+offset].set_levels(Zx, Wx)
        mesh.cs[ix+offset].bathy = B[ix]
        mesh.cs[ix+offset].update_geometry()
      
    if with_ghost:
        mesh.cs[2].crosssection_copy(mesh.cs[0])
        mesh.cs[2].crosssection_copy(mesh.cs[1])
        mesh.cs[X.size+1].crosssection_copy(mesh.cs[X.size+2])
        mesh.cs[X.size+1].crosssection_copy(mesh.cs[X.size+3])

    
def __setup_resampled_cross_sections__(mesh, X, B, Z, W, dx, dZmin, with_ghost):
    """ Setup cross-sections
    
        TODO : document arguments
    """
    
    # Compute resampled curvilinear abscissae
    L = X[-1] - X[0]
    nreg = int(math.ceil(L/dx))
    Lreg = dx * nreg
    Xreg = np.arange(X[0], X[0]+Lreg+dx, dx)
    #print(Xreg)
    
    # Resample variables
    nt = Z.shape[0]
    Breg = np.interp(Xreg, X, B)
    Zreg = np.zeros((nt, Xreg.size))
    Wreg = np.zeros((nt, Xreg.size))
    for it in range(0, nt):
        Zp = Z[it, :]
        Wp = W[it, :]
        valid_data = np.isfinite(Zp)
        Zp = Zp[valid_data]
        Xp = X[valid_data]
        spline = interpolate.InterpolatedUnivariateSpline(Xp, Zp, k=1)
        Zreg[it, :] = spline(Xreg)
        Wp = Wp[valid_data]
        Xp = X[valid_data]
        spline = interpolate.InterpolatedUnivariateSpline(Xp, Wp, k=1)
        Wreg[it, :] = spline(Xreg)
        indices = np.argwhere(Xreg < Xp[0]).flatten()
        if indices.size > 0:
            Wreg[it, indices] = np.maximum(Wreg[it, indices], Wp[0])
        indices = np.argwhere(Xreg > Xp[-1]).flatten()
        if indices.size > 0:
            Wreg[it, indices] = np.maximum(Wreg[it, indices], Wp[-1])
        
    __setup_cross_sections__(mesh, Xreg, Breg, Zreg, Wreg, dZmin, with_ghost)

    
def __finalize_mesh__(mesh, with_ghost):
    """ Setup cross-sections
    
        TODO : document arguments
    """
    
    slope_mean = 0.0
    if with_ghost:
        for ix in range(3, mesh.ncs-2):
        
            mesh.cs[ix].deltademi = mesh.cs[ix-1].x - mesh.cs[ix].x
            slope = (mesh.cs[ix-1].bathy - mesh.cs[ix].bathy) / abs(mesh.cs[ix].x - mesh.cs[ix-1].x)
            mesh.cs[ix-1].slope = slope
            slope_mean += slope
    else:
        for ix in range(1, mesh.ncs):
        
            mesh.cs[ix].deltademi = mesh.cs[ix-1].x - mesh.cs[ix].x
            slope = (mesh.cs[ix-1].bathy - mesh.cs[ix].bathy) / abs(mesh.cs[ix].x - mesh.cs[ix-1].x)
            mesh.cs[ix-1].slope = slope
            slope_mean += slope
        
    if with_ghost:
        slope_mean /= (mesh.ncs-5)
        mesh.cs[2].deltademi = mesh.cs[3].deltademi
        mesh.cs[0].slope = slope_mean
        mesh.cs[0].deltademi = mesh.cs[2].deltademi
        mesh.cs[1].slope = slope_mean
        mesh.cs[1].deltademi = mesh.cs[2].deltademi
        mesh.cs[mesh.ncs-3].slope = slope_mean
        #mesh.cs[mesh.ncs-3].deltademi = mesh.cs[mesh.ncs-4].deltademi
        mesh.cs[mesh.ncs-2].slope = slope_mean
        mesh.cs[mesh.ncs-2].deltademi = mesh.cs[mesh.ncs-3].deltademi
        mesh.cs[mesh.ncs-1].slope = slope_mean
        mesh.cs[mesh.ncs-1].deltademi = mesh.cs[mesh.ncs-3].deltademi
    else:
        slope_mean /= (mesh.ncs-1)
        mesh.cs[0].deltademi = mesh.cs[1].deltademi
        mesh.cs[mesh.ncs-1].slope = slope_mean
