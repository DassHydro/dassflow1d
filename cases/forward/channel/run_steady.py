import numpy as np
import os

import dassflow1d
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
from dassflow1d.post.results import plot_BZQ


def run_steady():
    
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
    # Set Strickler type (power_law_h:K=alpha*h^beta)
    mesh.set_strickler_type("powerlaw_h")
    # Set uniform parameters for Strickler
    mesh.set_uniform_strickler_parameters([25.0, 0.0])
    
    #------------------------------------------------------------------------------------------------------------------
    # MODEL SETUP
    #------------------------------------------------------------------------------------------------------------------
    # Initialise model from mesh
    model = m_sw_mono.Model(mesh)
    # Set upstream boundary condition (id:inflow discharge)
    model.bc[0].id = "discharge"
    # Set inflow timeseries
    model.bc[0].set_timeseries(t=[0, 1800, 3600, 5400, 7200], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    # Set downstream boundary condition (id:normal_depth)
    model.bc[1].id = "normal_depth"
    # Set scheme
    model.set_scheme("preissmann")
    
    #------------------------------------------------------------------------------------------------------------------
    # SIMULATION SETTINGS
    #------------------------------------------------------------------------------------------------------------------
    # Start time
    model.ts = 0
    # End time
    model.te = 0
    # Computation timestep
    model.dt = 60
    # Results timestep
    model.dtout = 600
    # Minimum water depth
    model.heps = 0.001
    # Implicit coefficient for Preissmann scheme
    model.theta_preissmann = 0.667

    #------------------------------------------------------------------------------------------------------------------
    # STEADY RUN
    #------------------------------------------------------------------------------------------------------------------
    print("=" * 80)
    print("Run steady simulation")
    print("=" * 80)
    # Run time loop
    model.steady_state()
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
    # Retrieve results field 'z' (water surface elevation) for segment 0 at output index 0
    z = model.get_segment_results(iseg=0, iout=0, field="z")
    # Retrieve results field 'q' (discharge) for segment 0 at output index 0
    q = model.get_segment_results(iseg=0, iout=0, field="q")
    # Create and show plot
    plot_BZQ(x, z, q, time_label="steady", bathy=bathy)


if __name__ == "__main__":
    run_steady()
