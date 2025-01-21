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
    mesh.set_uniform_strickler_parameters([20.0, 0.0])
    
    #------------------------------------------------------------------------------------------------------------------
    # MODEL SETUP
    #------------------------------------------------------------------------------------------------------------------
    # Initialise model from mesh
    model = m_sw_mono.Model(mesh)
    # Set upstream boundary condition (id:inflow discharge)
    model.bc[0].id = "discharge"
    # Set inflow timeseries
    model.bc[0].set_timeseries(t=[0, 1800, 3600, 5400, 7200], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    ## Set downstream boundary condition (id:normal_depth)
    model.bc[1].id = "normal_depth"
    # Set scheme
    model.set_scheme("preissmann")
    
    #------------------------------------------------------------------------------------------------------------------
    # SIMULATION SETTINGS
    #------------------------------------------------------------------------------------------------------------------
    # Start time
    model.ts = 0
    # End time
    model.te = 7200
    # Computation timestep
    model.dt = 60
    # Results timestep
    model.dtout = 600
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
    # VALIDATION
    #------------------------------------------------------------------------------------------------------------------
    print("=" * 80)
    print("Validation")
    print("=" * 80)
    # Retrieve field 'x' (distance from outlet) for segment 0
    x = mesh.get_segment_field(iseg=0, field="x")
    # Retrieve field 'bathy' (bathymetry) for segment 0
    bathy = mesh.get_segment_field(iseg=0, field="bathy")

    # Loop on output index to compute errors in L2 and Linf norms
    errZ_L2 = 0.0
    errZ_Linf = 0.0
    errQ_L2 = 0.0
    errQ_Linf = 0.0
    for iout in range(0, model.res.t.size):
            
        progress = iout * 100.0 / (model.res.t.size - 1)
        print("Validate results at time t=%f s ( %5.1f %% )" % (model.res.t[iout], progress))
        
        # Load DassFlow-1D v1 results
        if iout == 0:
            _, _, bathy_v1, _, q_v1, h_v1, _ = np.loadtxt("validation_data/result_initial.dat", unpack=True)
        elif iout == model.res.t.size-1:
            _, _, bathy_v1, _, q_v1, h_v1, _ = np.loadtxt("validation_data/result_final.dat", unpack=True)
        else:
            t = model.res.t[iout]
            _, _, bathy_v1, _, q_v1, h_v1, _ = np.loadtxt("validation_data/result_%12.6E.dat" % t, unpack=True)
        z_v1 = bathy_v1 + h_v1
            
        # Retrieve results field 'z' (water surface elevation) for segment 0
        z_v2 = model.get_segment_results(0, iout, "z")
        # Retrieve results field 'q' (discharge) for segment 0
        q_v2 = model.get_segment_results(0, iout, "q")
        
        # Make validation plot
        fig, axes = plt.subplots(2,1, sharex=True)
        axes[0].plot(x, bathy, 'k--', label="bathy")
        axes[0].plot(x, z_v2, 'b-', label="v2")
        axes[0].plot(x, z_v1, 'r+', label="v1")
        axes[0].legend()
        axes[1].plot(x, q_v2, 'b-', label="v2")
        axes[1].plot(x, q_v1, 'r+', label="v1")
        axes[1].legend()
        
        # Save plot
        plt.savefig(os.path.join("plot", "validation_%04i.png" % iout))
        plt.close(fig)
        
        # Print errors in L2 norm
        print("  - ErrZ_L2   = %12.6e, ErrQ_L2   = %12.6e" % (np.linalg.norm(z_v2-z_v1), np.linalg.norm(q_v2-q_v1)))
        print("  - ErrZ_Linf = %12.6e, ErrQ_Linf = %12.6e" % (np.max(np.abs(z_v2-z_v1)), np.max(np.abs(q_v2-q_v1))))
        
        # Append errors in L2 norm in lists and update errors in Linf norm
        errZ_L2 += np.linalg.norm(z_v2-z_v1)
        errZ_Linf = max(errZ_Linf, np.max(np.abs(z_v2-z_v1)))
        errQ_L2 += np.linalg.norm(z_v2-z_v1)
        errQ_Linf = max(errQ_Linf, np.max(np.abs(q_v2-q_v1)))

    print("Validation on whole case:")
    print("  - ErrZ_L2   = %12.6e, ErrQ_L2   = %12.6e" % (np.mean(errZ_L2), np.mean(errQ_L2)))
    print("  - ErrZ_LInf = %12.6e, ErrQ_LInf = %12.6e" % (errQ_Linf, errQ_Linf))
    

if __name__ == "__main__":
    run_unsteady()
