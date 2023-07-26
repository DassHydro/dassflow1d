import numpy as np
import os
import matplotlib.pyplot as plt

import dassflow1d
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
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
    # Set Strickler type (Power-law in h)
    mesh.set_strickler_type("powerlaw_h")
    # Set uniform parameters for Strickler
    mesh.set_uniform_strickler_parameters([40.0, 0.0])
    
    #------------------------------------------------------------------------------------------------------------------
    # MODEL SETUP
    #------------------------------------------------------------------------------------------------------------------
    # Initialise model from mesh
    model = m_sw_mono.Model(mesh)
    # Set upstream boundary condition (id:inflow discharge)
    model.bc[0].id = "discharge"
    # Load and set inflow timeseries
    model.bc[0].set_timeseries(t=[0.0, 3000.0], y=[100.0, 100.0])
    # Set downstream boundary condition (id:rating_curve)
    model.bc[1].id = "normal_depth"
    # Set scheme
    model.set_scheme("preissmann")
    
    #------------------------------------------------------------------------------------------------------------------
    # LATERAL INFLOW SETUP
    #------------------------------------------------------------------------------------------------------------------
    # Add inflow condition on segment 0 at coordinates (500.0, 0.0)
    # REMARK : when using DassFlow-1D v1.0 mesh files (old format, as for this case), the coordinates are the 2nd and 
    #          3rd values in the definition line of each cross-section
    model.add_inflow_condition(iseg=0, coords=[500.0, 0.0], t=[0, 3000.0], q=[20.0, 20.0])
    
    #------------------------------------------------------------------------------------------------------------------
    # SIMULATION SETTINGS
    #------------------------------------------------------------------------------------------------------------------
    # Start time
    model.ts = 0
    # End time
    model.te = 1500
    # Computation timestep
    model.dt = 5
    # Results timestep
    model.dtout = 50
    # Results output file
    model.output_file = "out/results.csv"
    # Minimum water depth
    model.heps = 0.001
    # Implicit coefficient for Preissmann scheme
    model.theta_preissmann = 0.667
    # Gravity constant
    model.gravity = 10.0

    #------------------------------------------------------------------------------------------------------------------
    # UNSTEADY RUN
    #------------------------------------------------------------------------------------------------------------------
    print("=" * 80)
    print("Run unsteady simulation")
    print("=" * 80)
    # Run time loop
    model.time_loop()
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
