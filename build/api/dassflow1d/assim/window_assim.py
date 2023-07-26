# encoding: utf-8
from __future__ import division, print_function, unicode_literals
import os
import matplotlib
if not "DISPLAY" in os.environ:
    matplotlib.use('Agg')
matplotlib.rc('text', usetex=True)
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

import dassflow1d
import dassflow1d.m_control as m_control
import dassflow1d.m_mesh as m_mesh

mincallback_dict = {}

def mincallback_window(x):

    global mincallback_dict
    
    # Increment iteration
    mincallback_dict["iteration"] +=1
    
    # Retrieve objects from dictionnary
    model = mincallback_dict["model"]
    control = mincallback_dict["control"]
    iteration = mincallback_dict["iteration"]
    x_prior = mincallback_dict["x_prior"]
    x_target = mincallback_dict["x_target"]
    xs_mesh = mincallback_dict["xs_mesh"]
    t_bc = mincallback_dict["t_bc"]
    
    # Retrieve real control vector
    x_current = np.zeros(x.size)
    control.get_real_control(x_current)
    
    if not os.path.isdir("plot"):
        os.mkdir("plot")
    plotdir = os.path.join("plot", "min_window")
    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)
        
    # Create figure
    fig, axes = plt.subplots(1, control.nitems, sharex=True, figsize=(control.nitems*6, 4))
        
    for iitem in range(0, control.nitems):
        
        # Retrieve parameters and variable to plot current control item
        if control.items[item_index].id[0:2] == b"BC":
            
            bc_index = int(control.items[item_index].id[2:5])
            if mdl%bc(index_bc+1)%id == "discharge":
                
                item_label = r"$Q$ ($m^3.s^{-1}$)"
                xlabel = "time (days)"
                item_xdata = model.bc[bc_index].t / 86400.0

        elif control.items[item_index].id[0:4] == b"KPAR":
            
            xlabel = "xs (km)"
            if xs_mesh is not None:
                xlabel = "xs (km)"
                item_xdata = xs_mesh / 1000.0
            else:
                xlabel = "index"
                item_xdata = np.arange(control.get_item_slice(iitem))
            par_index = int(control.items[item_index].id[4:6])
            if par_index == 1:
                if model.msh.strickler_type_code == m_mesh.strickler_type_constant:
                    item_label = r"$K$ ($m^{1/3}.s^{-1}$)"
                elif model.msh.strickler_type_code == m_mesh.strickler_type_powerlaw_h:
                    item_label = r"$\alpha$ ($m^{1/3}.s^{-1}$)"
                elif model.msh.strickler_type_code == m_mesh.strickler_type_einstein:
                    item_label = r"$KLOB$ ($m^{1/3}.s^{-1}$)"
                else:
                    item_label = control.items[item_index].id
            elif par_index == 2:
                if model.msh.strickler_type_code == m_mesh.strickler_type_powerlaw_h:
                    item_label = r"$\beta$ (-)"
                elif model.msh.strickler_type_code == m_mesh.strickler_type_einstein:
                    item_label = r"$KMIN$ ($m^{1/3}.s^{-1}$)"
                else:
                    item_label = control.items[item_index].id
            elif par_index == 3:
                if model.msh.strickler_type_code == m_mesh.strickler_type_einstein:
                    item_label = r"$KROB$ ($m^{1/3}.s^{-1}$)"
                else:
                    item_label = control.items[item_index].id
                    
        # TODO other possible items (BATHY, etc.)

        else:
            
            xlabel = "index"
            item_xdata = np.arange(control.get_item_slice(iitem))
            item_label = control.items[item_index].id

        # Create plot for current control item
        axes[iitem].plot(item_xdata, x_target[control.get_item_slice(iitem)], "r.", label="target")
        axes[iitem].plot(item_xdata, x_prior[control.get_item_slice(iitem)], "k--", label="prior")
        axes[iitem].plot(item_xdata, x_current[control.get_item_slice(iitem)], "b-", label="infered")
        axes[iitem].set_xlabel(xlabel)
        axes[iitem].set_ylabel(item_label)
    
    
    if final:
        iteration_in_title = "window %i - final" % window
        outfile = "plot/min_window/window%04i_final.png" % window
    else:
        iteration_in_title = "window %i - iteration %i" % (window, iteration)
        outfile = "plot/min_window/window%04i_iter%04i.png"  % (window, iteration)

    fig.suptitle("%s" % iteration_in_title)
      
    # Adjust plot layout
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Show or save figure
    plt.savefig(outfile)
    plt.close(fig)
    
    
        

def plot_estimations_vs_observations(obs):
    global plotsubdir

    plt.scatter(obs.obs[0, :], obs.est[0, :])
    xymin = min(np.min(obs.obs[0, :]), np.min(obs.est[0, :]))
    xymax = max(np.max(obs.obs[0, :]), np.max(obs.est[0, :]))
    plt.plot([xymin, xymax], [xymin, xymax], 'k--')
    plt.ylim((xymin, xymax))
    plt.xlim((xymin, xymax))
    #plt.savefig("plot/%s/obs_vs_est.png" % plotsubdir)
    #plt.close(plt.gcf())
    plt.show()

def window_assim(model, control, obs, window, xtarget=None, minfunc=None, minargs=(), minkwargs={}):
    global mincallback_dict

    if minfunc is None:
        
        minfunc = scipy.optimize.minimize
        minargs = (dassflow1d.calc_cost_and_gradients_scipy_minimize, control.x)
        minkwargs = {"args" : (model, control, obs),
                     "jac" : True,
                     "tol" : 1e-3, 
                     "method" : 'L-BFGS-B',
                     "options" : {"disp": True},
                     "callback" : mincallback_window}
        
    # Store model start and end times
    ts = model.ts
    te = model.te
        
    # Set initial window index and start time
    window_index = 1
    tws = model.ts
    
    # Loop until window start time exceed model end time
    while tws < te:
        
        # Set window end time and update model start and end time with window times
        twe = min(tws+window, te)
        model.ts = tws
        model.te = twe
        
        window_obs = obs.extract_window(tws, twe)
        
        # Retrieve prior
        xprior = np.zeros(control.x.size)
        control.get_real_control(xprior)
        
        if model.msh.nseg == 1:
            xs_mesh = model.msh.get_segment_field(0, "x")
        else:
            # TODO multiple segments ?
            xs_mesh = None
        
        # Set callback dictionnary
        mincallback_dict = {"iteration" : 1,
                            "model" : model,
                            "control" : control,
                            "x_prior" : xprior,
                            "x_target": xtarget,
                            "xs_mesh": xs_mesh,
                            "t_bc": []}
                            
        window_index = 1
        tws = model.ts
        
        # Run minimization
        print("=" * 80)
        print(" MINIMIZATION on WINDOW %i (%.2f - %.2f)" % (window_index, model.ts / 86400.0, model.te / 86400.0))
        print("=" * 80)
        iteration = 0
        res = minfunc(*minargs, **minkwargs)
        
        
        plot_estimations_vs_observations(obs)
        
        print("")
        print("status=", model.status)
        iteration = -1
        mincallback_dict["iteration"] = -1
        mincallback_window(res.x)
        
        # Transfer infered control vector for next window
        control.set_prior(res.x)
        
        # Disable initialisation (last state is used as initialisation for the next window)
        model.external_initialisation = True
        
        # Increment window index
        window_index += 1
