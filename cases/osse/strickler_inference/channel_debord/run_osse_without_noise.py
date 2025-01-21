import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import scipy.optimize

import dassflow1d
import dassflow1d.m_control as m_control
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_obs as m_obs

control = None
x_target = None
x_prior = None
iteration = 0
ite_cost = []
ite_x = []


def plot_control(ite_x, x_target):

    global control
    
    # Create figure
    fig, axes = plt.subplots(3, 1, sharex=True)

    # Make subplots of KLOB (Left Overbank), KMC (Main Channel) and KROB (Right Overbank) w.r.t iterations 
    labels = ["KLOB", "KMC", "KROB"]
    for i in range(0, 3):
        axes[i].plot(np.arange(1, len(ite_x)+1), np.ones(len(ite_x)) * x_target[i], "r--", label="target")
        axes[i].plot(np.arange(1, len(ite_x)+1), ite_x, "b-", label="infered")
        axes[i].set_ylabel(labels[i])
        axes[i].xaxis.set_major_locator(MaxNLocator(integer=True))
    axes[2].set_xlabel("iterations")
    axes[2].legend()
        
    # Save plot to file
    plt.savefig("plot/x_without_noise.png")
    plt.close(plt.gcf())


def plot_cost(ite_cost):

    # Make plot of cost w.r.t iterations 
    plt.plot(np.arange(1, len(ite_cost)+1), ite_cost, "g-")
    plt.xlabel("iteration")
    plt.ylabel("cost (L2-norm)")
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        
    # Save plot to file
    plt.savefig("plot/cost_without_noise.png")
    plt.close(plt.gcf())


def callback_minimize(x, final=False):
    
    global control
    global iteration
    global ite_cost
    global ite_x
    global model
    global x_prior
    global x_target
    global xs_mesh
    
    # Increment iteration counter
    iteration+=1
    
    if final:
        
        # Plot control items
        plot_control(ite_x, x_target)

        # Plot cost
        plot_cost(ite_cost)
        
    else:
        
        # Store control vector and cost
        # Retrieve real control values (reverse the change of variable if activated)
        x_iter = np.zeros(control.x.size)
        control.get_real_control(x_iter)
        print("TARGET :", x_target)
        print("INFERED: ", x_iter)

        # Store real control values and cost
        ite_x.append(x_iter)
        ite_cost.append(model.__last_cost__)


def run_osse():
    
    #------------------------------------------------------------------------------------------------------------------
    # GLOBAL VALUES (to make them accessible for the callback_minimize function)
    #------------------------------------------------------------------------------------------------------------------
    global control
    global model
    global obs
    global x_target
    global x_prior
    global xs_mesh
    
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
    # Set strickler type (Debord formula)
    mesh.set_strickler_type("Debord")
    # Set uniform parameters for Strickler
    mesh.set_strickler_fields_segment([[10.0], [20.0], [10.0]])
    # Retrieve curvilinear abscissae for plots
    xs_mesh = mesh.get_segment_field(iseg=0, field="x")
    
    #------------------------------------------------------------------------------------------------------------------
    # MODEL SETUP
    #------------------------------------------------------------------------------------------------------------------
    # Initialise model from mesh
    model = m_sw_mono.Model(mesh)
    # Set upstream boundary condition (id:inflow discharge)
    model.bc[0].id = "discharge"
    # Set inflow timeseries
    t1 = np.linspace(0, 864000, 21, endpoint=True)
    qmin=60
    qmax=80
    q1 = qmax - 2 * (qmax-qmin) * abs(t1 - 432000.0) / 864000.0
    model.bc[0].set_timeseries(t=t1, y=q1)
    # Set downstream boundary condition (id:normal_depth)
    model.bc[1].id = "normal_depth"
    # Set scheme
    model.set_scheme("preissmann")
    
    #------------------------------------------------------------------------------------------------------------------
    # SIMULATION SETTINGS
    #------------------------------------------------------------------------------------------------------------------
    # Set start time
    model.ts = 0
    # Set end time
    model.te = 864000
    # Set computation timestep
    model.dt = 300
    # Set output timestep to -1 (no output)
    model.dtout = -1
    # Set minimal depth
    model.heps = 0.01
    # Set gamma_reg to zero (no regularisation)
    model.gamma_reg = 0.0
    # Set implicit theta for Preissmann scheme
    model.theta_preissmann = 0.667

    #------------------------------------------------------------------------------------------------------------------
    # GENERATION OF OBSERVATIONS
    #------------------------------------------------------------------------------------------------------------------
    print("=" * 80)
    print(" GENERATION OF OBSERVATIONS")
    print("=" * 80)
    # Set observations at every cross-section at regular timestep
    obs = m_obs.all_observed(np.linspace(10800.0, 864000, 800, endpoint=True), mesh)
    # Run direct model to generate observations
    model.generate_osse_observations(obs)
    print("")

    #------------------------------------------------------------------------------------------------------------------
    # CONTROL VECTOR INITIALISATION
    #------------------------------------------------------------------------------------------------------------------
    # Init control for minimization
    control = m_control.Control()
    # Add strickler parameters in control
    control.add_strickler_in_control(model, mesh)
    # Copy current control (values used for the generation of observations) to x_target (for plots)
    x_target = control.x.copy()
    # Set prior values
    control.x0[control.get_item_slice(0)] = 15.0
    control.x0[control.get_item_slice(1)] = 25.0
    control.x0[control.get_item_slice(2)] = 15.0
    # Copy prior control to x_prior (for plots)
    x_prior = control.x0.copy()
    # Copy prior control to current control
    control.x[:] = control.x0[:]
    
    #------------------------------------------------------------------------------------------------------------------
    # MINIMIZATION
    #------------------------------------------------------------------------------------------------------------------
    print("=" * 80)
    print(" MINIMIZATION")
    print("=" * 80)
    # Run minimisation
    res = scipy.optimize.minimize(dassflow1d.calc_cost_and_gradients_scipy_minimize, control.x, 
                                  args=(model, control, obs), jac=True, tol = 1e-7, 
                                  method='L-BFGS-B', options={"disp": True},
                                  callback=callback_minimize)
    # Use minimisation callback with iteration=-1 to make final plots
    iteration = -1
    callback_minimize(res.x, final=True)
    print("")


if __name__ == "__main__":
    run_osse()
