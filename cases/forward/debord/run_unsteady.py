import numpy as np
import os
import matplotlib.pyplot as plt

import dassflow1d
#import dassflow1d.m_control as m_control
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
#import dassflow1d.m_obs as m_obs
from dassflow1d.post.results import plot_BZQ


def run_unsteady():
    
    #------------------------------------------------------------------------------------------------------------------
    # INIT
    #------------------------------------------------------------------------------------------------------------------
    # Create non-existing directories
    if not os.path.isdir("out"):
        os.mkdir("out")
    if not os.path.isdir("plot"):
        os.mkdir("plot")

    #------------------------------------------------------------------------------------------------------------------
    # MESH SETUP
    #------------------------------------------------------------------------------------------------------------------
    # Read mesh
    mesh = dassflow1d.read_mesh("mesh.geo")
    # Set Strickler type (Debord formula)
    mesh.set_strickler_type("Debord")
    # Read Strickler values from file
    mesh.read_spatial_field("ks.dat")
    
    #------------------------------------------------------------------------------------------------------------------
    # MODEL SETUP
    #------------------------------------------------------------------------------------------------------------------
    # Initialise model from mesh
    model = m_sw_mono.Model(mesh)
    # Set upstream boundary condition (id:inflow discharge)
    model.bc[0].id = "discharge"
    # Load and set inflow timeseries
    t, q = np.loadtxt("upstream.csv", delimiter=",", skiprows=1, unpack=True)
    model.bc[0].set_timeseries(t=t, y=q)
    # Set downstream boundary condition (id:rating_curve)
    model.bc[1].id = "rating_curve"
    # Load and set rating curve table
    z, q = np.loadtxt("downstream.csv", delimiter=",", skiprows=1, unpack=True)
    model.bc[1].set_rating_curve(z=z, q=q)
    # Set scheme
    model.set_scheme("preissmann")
    
    #------------------------------------------------------------------------------------------------------------------
    # SIMULATION SETTINGS
    #------------------------------------------------------------------------------------------------------------------
    # Start time
    model.ts = 0
    # End time
    model.te = 43200
    # Computation timestep
    model.dt = 300
    # Results timestep
    model.dtout = 3600
    # Results output file
    model.output_file = "out/results.csv"
    # Minimum water depth
    model.heps = 0.001
    # Implicit coefficient for Preissmann scheme
    model.theta_preissmann = 0.667

    #------------------------------------------------------------------------------------------------------------------
    # UNSTEADY RUN
    #------------------------------------------------------------------------------------------------------------------
    print("=" * 80)
    print("Run unsteady simulation")
    print("=" * 80)
    model.run_unsteady()
    print("")

    #------------------------------------------------------------------------------------------------------------------
    # POST-PROCESSING
    #------------------------------------------------------------------------------------------------------------------
    print("=" * 80)
    print("Post-processing")
    print("=" * 80)
    # Retrieve field 'x' (distance from outlet) for segment 0
    x = mesh.get_segment_field(iseg=0, field="x")
    # Retrieve field 'bathy' (bathymetry) for segment 0
    bathy = mesh.get_segment_field(iseg=0, field="bathy")
    # Loop on output index
    for iout in range(0, len(model.res.t)):
        
        progress = iout * 100.0 / (model.res.t.size - 1)
        print("Plot results at time t=%f s ( %5.1f %% )" % (model.res.t[iout], progress))
        
        # Retrieve time for current output index
        t = model.get_results_time(iout=iout)
        # Retrieve results field 'q' (discharge) for segment 0 and current output index
        q = model.get_segment_results(iseg=0, iout=iout, field="q")
        # Retrieve results field 'z' (water surface elevation) for segment 0 and current output index
        z = model.get_segment_results(iseg=0, iout=iout, field="z")
            
        # Make plot
        plot_BZQ(x, z, q, "t=%.1f s" % t, bathy=bathy, outfile="plot/results_%04i.png" % iout)
            

if __name__ == "__main__":
    
    run_unsteady()
