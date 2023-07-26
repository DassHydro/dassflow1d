# encoding: utf-8
from __future__ import division, print_function, unicode_literals
import os
import matplotlib
if not "DISPLAY" in os.environ:
    matplotlib.use('Agg')
matplotlib.rc('text', usetex=True)
import numpy as np
import matplotlib.pyplot as plt


def plot_BZQ(xs, z, q, time_label, bathy=None, xs_unit=None, river_name=None, outfile=None):
  
    # CHECK-UP and type cast
    if not isinstance(xs, np.ndarray):
        xs = np.array(xs)
    if not isinstance(z, np.ndarray):
        z = np.array(z)
    if not isinstance(q, np.ndarray):
        q = np.array(q)
    if not isinstance(time_label, str):
      raise ValueError("'time_label' must be a string")
    if not isinstance(bathy, np.ndarray):
        bathy = np.array(bathy)
    if river_name is not None:
        if not isinstance(river_name, str):
            raise ValueError("'river_name' must be a string")
    if outfile is not None:
        if not isinstance(outfile, str):
            raise ValueError("'outfile' must be a string")
    
    # Compute units and multiplier for curvilinear abscissae
    if xs_unit is None:
        if np.max(xs) < 2000:
            xs_unit = "m"
        else:
            xs_unit = "km"
    if xs_unit == "m":
        xs_mult = 1.0
    else:
        xs_mult = 0.001
      
    # Initialise figure and subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    
    # Plot bathymetry
    if bathy is not None:
        if isinstance(xs, list):
            for i in range(0, len(xs)):
                ax1.plot(xs[i]*xs_mult, bathy[i], 'k-', alpha=0.5)
        else:
            ax1.fill_between(xs*xs_mult, np.min(bathy), bathy, facecolor='k', alpha=0.5)
        
    # Plot free surface height
    if isinstance(xs, list):
        for i in range(0, len(xs)):
            ax1.plot(xs[i]*xs_mult, z[i], 'b-')
    else:
        ax1.plot(xs*xs_mult, z, 'b-')
    
    # Set axes labels
    ax1.set_xlabel(r"$x$ $(%s)$" % xs_unit)
    ax1.set_ylabel(r"$H$ $(m)$")
    
    # Plot discharge
    if isinstance(xs, list):
        for i in range(0, len(xs)):
            ax1.plot(xs[i]*xs_mult, q[i], 'r-')
    else:
        ax2.plot(xs*xs_mult, q, 'r-')
    
    # Set axes labels
    ax2.set_xlabel(r"$x$ $(%s)$" % xs_unit)
    ax2.set_ylabel(r"$Q$ $(m^3.s^{-1})$")
    
    # Set main title
    if river_name is not None:
      fig.suptitle("%s - %s" % (river_name, time_label))
    else:
      fig.suptitle("%s" % time_label)
      
    # Adjust plot layout
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Show or save figure
    if outfile:
        plt.savefig(outfile)
        plt.close(fig)
    else:
        plt.show()
        
        
def plot_hQ(xs, h, q, time_label, xs_unit=None, river_name=None, outfile=None):
  
    # CHECK-UP and type cast
    if not isinstance(xs, np.ndarray):
        xs = np.array(xs)
    if not isinstance(h, np.ndarray):
        h = np.array(z)
    if not isinstance(q, np.ndarray):
        q = np.array(q)
    if not isinstance(time_label, str):
      raise ValueError("'time_label' must be a string")
    if river_name is not None:
        if not isinstance(river_name, str):
            raise ValueError("'river_name' must be a string")
    if outfile is not None:
        if not isinstance(outfile, str):
            raise ValueError("'outfile' must be a string")
    
    # Compute units and multiplier for curvilinear abscissae
    if xs_unit is None:
        if np.max(xs) < 2000:
            xs_unit = "m"
        else:
            xs_unit = "km"
    if xs_unit == "m":
        xs_mult = 1.0
    else:
        xs_mult = 0.001
      
    # Initialise figure and subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    
    # Plot water depth
    ax1.plot(xs*xs_mult, h, 'b-')
    
    # Set axes labels
    ax1.set_xlabel(r"$x$ $(%s)$" % xs_unit)
    ax1.set_ylabel(r"$h$ $(m)$")
    
    # Plot discharge
    ax2.plot(xs*xs_mult, q, 'r-')
    
    # Set axes labels
    ax2.set_xlabel(r"$x$ $(%s)$" % xs_unit)
    ax2.set_ylabel(r"$Q$ $(m^3.s^{-1})$")
    
    # Set main title
    if river_name is not None:
      fig.suptitle("%s - %s" % (river_name, time_label))
    else:
      fig.suptitle("%s" % time_label)
      
    # Adjust plot layout
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Show or save figure
    if outfile:
        plt.savefig(outfile)
        plt.close(fig)
    else:
        plt.show()


def plot_BZQW(xs, z, q, w, time_label, bathy=None, a=None, xs_unit=None, river_name=None, outfile=None):
  
    # CHECK-UP and type cast
    if not isinstance(xs, np.ndarray):
        xs = np.array(xs)
    if not isinstance(z, np.ndarray):
        z = np.array(z)
    if not isinstance(q, np.ndarray):
        q = np.array(q)
    if not isinstance(w, np.ndarray):
        w = np.array(w)
    if not isinstance(time_label, str):
      raise ValueError("'time_label' must be a string")
    if not isinstance(bathy, np.ndarray):
        bathy = np.array(bathy)
    if river_name is not None:
        if not isinstance(river_name, str):
            raise ValueError("'river_name' must be a string")
    if outfile is not None:
        if not isinstance(outfile, str):
            raise ValueError("'outfile' must be a string")
    
    # Compute units and multiplier for curvilinear abscissae
    if xs_unit is None:
        if np.max(xs) < 2000:
            xs_unit = "m"
        else:
            xs_unit = "km"
    if xs_unit == "m":
        xs_mult = 1.0
    else:
        xs_mult = 0.001
      
    # Initialise figure and subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    
    # Plot bathymetry
    if bathy is not None:
        ax1.fill_between(xs*xs_mult, np.min(bathy), bathy, facecolor='k', alpha=0.5)
        
    # Plot free surface height
    ax1.plot(xs*xs_mult, z, 'b-')
    
    # Set axes labels
    #ax1.set_xlabel(r"$x$ $(%s)$" % xs_unit)
    ax1.set_ylabel(r"$H$ $(m)$")
    
    # Plot discharge
    ax2.plot(xs*xs_mult, q, 'r-')
    if a is not None:
        ax2twin = ax2.twinx()
        ax2twin.plot(xs*xs_mult, a, 'g-')
    
    # Set axes labels
    #ax2.set_xlabel(r"$x$ $(%s)$" % xs_unit)
    ax2.set_ylabel(r"$Q$ $(m^3.s^{-1})$")
    
    # Plot wave speeds
    ax3.plot(xs*xs_mult, w, 'g-')
    
    # Set axes labels
    ax3.set_xlabel(r"$x$ $(%s)$" % xs_unit)
    ax3.set_ylabel(r"$w$ $(m)$")
    
    # Set main title
    if river_name is not None:
      fig.suptitle("%s - %s" % (river_name, time_label))
    else:
      fig.suptitle("%s" % time_label)
      
    # Adjust plot layout
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Show or save figure
    if outfile:
        plt.savefig(outfile)
        plt.close(fig)
    else:
        plt.show()


def plot_BZQwaves(xs, z, q, u, c, time_label, bathy=None, a=None, xs_unit=None, river_name=None, outfile=None):
  
    # CHECK-UP and type cast
    if not isinstance(xs, np.ndarray):
        xs = np.array(xs)
    if not isinstance(z, np.ndarray):
        z = np.array(z)
    if not isinstance(q, np.ndarray):
        q = np.array(q)
    if not isinstance(u, np.ndarray):
        u = np.array(u)
    if not isinstance(c, np.ndarray):
        c = np.array(c)
    if c.size == 1:
        c = np.ones(u.size) * c
    if not isinstance(time_label, str):
      raise ValueError("'time_label' must be a string")
    if not isinstance(bathy, np.ndarray):
        bathy = np.array(bathy)
    if river_name is not None:
        if not isinstance(river_name, str):
            raise ValueError("'river_name' must be a string")
    if outfile is not None:
        if not isinstance(outfile, str):
            raise ValueError("'outfile' must be a string")
    
    # Compute units and multiplier for curvilinear abscissae
    if xs_unit is None:
        if np.max(xs) < 2000:
            xs_unit = "m"
        else:
            xs_unit = "km"
    if xs_unit == "m":
        xs_mult = 1.0
    else:
        xs_mult = 0.001
      
    # Initialise figure and subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    
    # Plot bathymetry
    if bathy is not None:
        ax1.fill_between(xs*xs_mult, np.min(bathy), bathy, facecolor='k', alpha=0.5)
        
    # Plot free surface height
    ax1.plot(xs*xs_mult, z, 'b-')
    
    # Set axes labels
    #ax1.set_xlabel(r"$x$ $(%s)$" % xs_unit)
    ax1.set_ylabel(r"$H$ $(m)$")
    
    # Plot discharge
    ax2.plot(xs*xs_mult, q, 'r-')
    if a is not None:
        ax2twin = ax2.twinx()
        ax2twin.plot(xs*xs_mult, a, 'g-')
    
    # Set axes labels
    #ax2.set_xlabel(r"$x$ $(%s)$" % xs_unit)
    ax2.set_ylabel(r"$Q$ $(m^3.s^{-1})$")
    
    # Plot wave speeds
    ax3.plot(xs*xs_mult, u+c, 'b-')
    ax3.plot(xs*xs_mult, u, 'r--')
    ax3.plot(xs*xs_mult, c, 'g--')
    
    # Set axes labels
    ax3.set_xlabel(r"$x$ $(%s)$" % xs_unit)
    ax3.set_ylabel(r"$u,c,w$ $(m.s^{-1})$")
    
    # Set main title
    if river_name is not None:
      fig.suptitle("%s - %s" % (river_name, time_label))
    else:
      fig.suptitle("%s" % time_label)
      
    # Adjust plot layout
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    
    # Show or save figure
    if outfile:
        plt.savefig(outfile)
        plt.close(fig)
    else:
        plt.show()
