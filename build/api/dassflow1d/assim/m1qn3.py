# encoding: utf-8
from __future__ import division, print_function, unicode_literals
import numpy as np

from dassflow1d.m_minimization import M1Qn3Comm, m1qn3_reverse_mode
from scipy.optimize import OptimizeResult


def minimize_m1qn3(fun, x, jac=False, status=False, tol=1e-6, callback=None, **kwargs):
    
    if jac == False:
        raise NotImplementedError("M1QN3 with 'jac' is not True is not implemented")
    
    
    # Setup parameters
    n = len(x)
    args = ()
    if "args" in kwargs:
        args = kwargs["args"]
    normtype = kwargs["normtype"] if "normtype" in kwargs else 'two'
    
    comm = M1Qn3Comm(n)
    comm.dxmin = tol
    if "disp" in kwargs:
        if disp is True:
            comm.impres = 3
        elif disp is False:
            comm.impres = 0
        else:
            raise ValueError("'disp' must be True or False")
    if "epsg" in kwargs:
        comm.epsg = kwargs["epsg"]
    if "nsim" in kwargs:
        comm.nsim = kwargs["nsim"]
    comm.reverse = 1
    comm.indic = 4
    #idz = np.zeros(5, dtype=int)
    #dz = np.zeros(ndz, dtype=int)
        
    while comm.reverse >= 0 and comm.indic != 0:
        
        if comm.indic == 4:
            
            # M1QN3 requires new estimates of the cost and gradients
            if status:
                f, g, sim_status = fun(x, *args)
            else:
                f, g = fun(x, *args)
                sim_status = 0
            
        if sim_status != 0:
            comm.indic = -1
            
        comm.df1 = 0.5 * f
    
        m1qn3_reverse_mode(n, x, f, g, comm)
        
        
        if comm.indic == 1 and callback:
            callback(x)
            
    res = OptimizeResult.fromkeys(["x", "success", "status", "fun", "jac", "nfev", "nit"])
    res["x"] = x
    res["success"] = comm.omode == 1
    res["status"] = comm.omode
    res["fun"] = f
    res["jac"] = g
    res["nfev"] = comm.nsim
    res["nit"] = comm.niter

    return res
