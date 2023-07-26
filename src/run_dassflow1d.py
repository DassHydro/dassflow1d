import argparse
import json
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize

import dassflow1d
import dassflow1d.m_control as m_control
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_obs as m_obs


def load_configuration(fname):
    
    if os.path.splitext(fname)[1] == ".json":
        
        # Load configuration from json file
        with open(fname) as fin:
            config = json.load(fin)
            
    else:
        
        raise RuntimeError("configuration file must be a json file")
    
    
    # TODO check config content
    
    return config


def load_timeseries(fname, datestart=None):
    
    if os.path.splitext(fname)[1] == ".csv":
        
        if datestart is None:
            
            data = pd.read_csv(fname, sep=";")
            return data.iloc[:, 0].values, data.iloc[:, 1].values
        
        else:
            
            data = pd.read_csv(fname, sep=";", parse_dates=[0])
            t = ((data.iloc[:, 0].dt.tz_localize(None) - np.datetime64(datestart)) / np.timedelta64(1, "s")).values
            
            return t, data.iloc[:, 1].values
            
    else:
        
        raise RuntimeError("Timeseries file must be a CSV file")


def create_model(config):
    """ Setup model
    
        Parameters
        ----------
            config: dict
                Run configuration
                
        Return
        ------
            m_sw_mono.Model
                DassFlow-1D model
    """

    print("=" * 80)
    print(" SETUP MODEL")
    print("=" * 80)
    
    # Load mesh
    mesh_fname = config["model"]["mesh_file"]
    mesh = dassflow1d.read_mesh(mesh_fname)
    
    # Resample mesh if requested
    if "mesh_resample" in config["model"]:
        mesh.resample(float(config["model"]["mesh_resample"]))
        
    # Configure Strickler type
    if "strickler_type" in config["model"]:
        mesh.set_strickler_type(config["model"]["strickler_type"])
        
    # Setup strickler values
    if config["model"]["strickler_values"]["type"] == "uniform":
        
        # Set uniform values
        values = config["model"]["strickler_values"]["values"]
        mesh.set_uniform_strickler_parameters(values)
        
    elif config["model"]["strickler_values"]["type"] == "fields":
        pass
        #if "strickler_fields" in config["model"]:
            #if config["model"]["strickler_fields"] == "segment":
                #pass
        #else:
            #raise RuntimeError("'strickler_fields' parameter not specified in 'model' section")
    else:
        raise RuntimeError("wrong type of strickler values : %s" % config["model"]["strickler_values"]["type"])
    
    # Create model from mesh
    model = m_sw_mono.Model(mesh)
    
    # Retrieve simulation dates or duration (in seconds)
    if "date_start" in config["model"]:
        date_start = config["model"]["date_start"]
    else:
        date_start = None
    if "date_end" in config["model"]:
        date_end = config["model"]["date_end"]
    else:
        date_end = None
    if "duration" in config["model"]:
        duration = config["model"]["duration"]
    else:
        duration = None
    
    # Setup boundary conditions
    for ibc, bc_config in enumerate(config["model"]["boundary_conditions"]):
        model.bc[ibc].id = bc_config["id"]
        if "timeseries_file" in bc_config:
            t, y = load_timeseries(bc_config["timeseries_file"], datestart=date_start)
            model.bc[ibc].set_timeseries(t, y)
    
    # Setup inflow conditions
    if "inflow_conditions" in config["model"]:
        for iic, ic_config in enumerate(config["model"]["inflow_conditions"]):
            iseg = ic_config["segment"]
            if "coords" in ic_config:
                coords = ic_config["coords"]
            else:
                raise RuntimeError("'coords' must be specified for inflow conditions")
            if "timeseries_file" in ic_config:
                t, y = load_timeseries(ic_config["timeseries_file"], datestart=date_start)
            else:
                raise RuntimeError("'timeseries_file' must be specified for inflow conditions")
            model.add_inflow_condition(iseg=iseg, coords=coords, t=t, q=y)
            
    # Set model scheme
    if "numerical_scheme" in config["model"]:
        model.set_scheme(config["model"]["numerical_scheme"])
    else:
        model.set_scheme("preissmann")

    # Set start time
    model.ts = 0
    
    # Set end time
    if duration is None:
        model.te = (np.datetime64(date_end) - np.datetime64(date_start)) / np.timedelta64(1, "s")
    else:
        model.te = duration
        
    # Set computation timestep
    model.dt = config["model"]["timestep"]
    
    # Set output timestep
    model.dtout = config["model"]["output_timestep"]
    model.output_resampled_cs = True
    
    # Set output file
    if "output_file" in config["model"]:
        model.output_file = config["model"]["output_file"]
    else:
        if not os.path.isdir("out"):
            os.mkdir("out")
        model.output_file = "out/results.csv"
    
    # Set minimal depth
    if "heps" in config["model"]:
        model.heps = config["model"]["heps"]
    else:
        model.heps = 0.01
    
    # Set minimal discharge
    if "qeps" in config["model"]:
        model.qeps = config["model"]["qeps"]
    else:
        model.qeps = 0.1
        
    ## Set gamma for the regularization
    #model.gamma_reg = 0.0
    
    # Set implicit theta for Preissmann scheme
    if "theta_preissmann" in config["model"]:
        model.theta_preissmann = config["model"]["theta_preissmann"]
    else:
        model.theta_preissmann = 0.667
        
            
    return model


def load_observations(config, model):
    """ Load observations
    
        Parameters
        ----------
            config: dict
                Run configuration
            model: m_sw_mono.Model
                DassFlow-1D model
                
        Return
        ------
            m_obs.Observations
                Observations object
    """

    print("=" * 80)
    print(" SETUP OBSERVATIONS STATIONS")
    print("=" * 80)
    
    # Retrieve date start if specified
    if "date_start" in config["model"]:
        date_start = config["model"]["date_start"]
    else:
        date_start = None
    
    # Allocate observations object
    nobs = len(config["observations"])
    obs = m_obs.Observations(nobs)
    
    # Setup Verdun-sur-Garonne station
    for iobs, obs_config in enumerate(config["observations"]):
        
        # Load timeseries
        tobs, Hobs = load_timeseries(obs_config["timeseries_file"], datestart=date_start)
        
        # Create array of observed H and W, restricted to the simulation window
        Hobs = Hobs[tobs <= model.te]
        tobs = tobs[tobs <= model.te]
        HWobs = np.ones((2, Hobs.size)) * -1e+99
        HWobs[0, :] = Hobs
        
        # Setup station
        obs.stations[iobs].setup(model.msh, tobs, HWobs, indices=obs_config["index"])
        
    print("")
    
    return obs


def load_control(config, model):
    """ Setup control
    
        Parameters
        ----------
            config: dict
                Run configuration
            model: m_sw_mono.Model
                DassFlow-1D model
                
        Return
        ------
            m_control.Control
                Control object
    """

    print("=" * 80)
    print(" SETUP CONTROL")
    print("=" * 80)

    # Create control object
    control = m_control.Control()
    
    # Retrieve model mesh
    mesh = model.msh
    
    # Add items in control
    for item_config in config["control"]:
        
        if item_config["variable"] == "K":

            # Check that Strickler type is constant
            if model.msh.strickler_type_code != m_mesh.strickler_type_constant:
                raise RuntimeError("Cannot add variable 'K' in control if strickler type is not 'constant'")
            
            # Add Strickler in control
            control.add_strickler_component_in_control(model, mesh, 0)

        elif item_config["variable"] == "alpha":
            
            # Check that Strickler type is powerlaw_h
            if model.msh.strickler_type_code != m_mesh.strickler_type_powerlaw_h:
                raise RuntimeError("Cannot add variable 'alpha' in control if strickler type is not 'powerlaw_h'")
            
            # Add alpha
            control.add_strickler_component_in_control(model, mesh, 0)

        elif item_config["variable"] == "beta":
            
            # Check that Strickler type is powerlaw_h
            if model.msh.strickler_type_code != m_mesh.strickler_type_powerlaw_h:
                raise RuntimeError("Cannot add variable 'beta' in control if strickler type is not 'powerlaw_h'")
            
            # Add beta
            control.add_strickler_component_in_control(model, mesh, 1)

        elif item_config["variable"] == "KLOB":
            
            # Check that Strickler type is Debord
            if model.msh.strickler_type_code != m_mesh.strickler_type_debord:
                raise RuntimeError("Cannot add variable 'KLOB' in control if strickler type is not 'Debord'")
            
            # Add KLOB
            control.add_strickler_component_in_control(model, mesh, 0)

        elif item_config["variable"] == "KCH":
            
            # Check that Strickler type is Debord
            if model.msh.strickler_type_code != m_mesh.strickler_type_debord:
                raise RuntimeError("Cannot add variable 'KCH' in control if strickler type is not 'Debord'")
            
            # Add KCH
            control.add_strickler_component_in_control(model, mesh, 1)

        elif item_config["variable"] == "KROB":
            
            # Check that Strickler type is Debord
            if model.msh.strickler_type_code != m_mesh.strickler_type_debord:
                raise RuntimeError("Cannot add variable 'KROB' in control if strickler type is not 'Debord'")
            
            # Add KROB
            control.add_strickler_component_in_control(model, mesh, 2)
            
        else:
            
            raise RuntimeError("Cannot add unknown variable '%s' in control" % item_config["variable"])

    print("Control size: %i" % control.x.size)
    print("")
    
    return control


def run_direct(config, display_internal_counters=False):
  
    # Create model
    model = create_model(config)
    
    # Run direct
    print("=" * 80)
    print(" SIMULATION")
    print("=" * 80)
    
    # Run direct model
    model.time_loop()
    if display_internal_counters:
        print("- Number of overbank flow occurences: %i" % model.internal_counters[0])
        print("  - Left overbank flow occurences : %i" % model.internal_counters[1])
        print("  - Right overbank flow occurences: %i" % model.internal_counters[2])
        print("- Number of overflow (flow above highest elevation) occurences: %i" % model.internal_counters[3])
    print("")


def run_calc_cost(config, display_internal_counters=False):
  
    # Create model
    model = create_model(config)
  
    # Load observations
    obs = load_observations(config, model)
  
    # Setup dummy control
    control = m_control.Control()
    
    # Run direct 
    print("=" * 80)
    print(" SIMULATION")
    print("=" * 80)
    
    # Run direct model to compute cost
    cost = dassflow1d.calc_cost(model, control, obs)
    if display_internal_counters:
        print("- Number of overbank flow occurences: %i" % model.internal_counters[0])
        print("  - Left overbank flow occurences : %i" % model.internal_counters[1])
        print("  - Right overbank flow occurences: %i" % model.internal_counters[2])
        print("- Number of overflow (flow above highest elevation) occurences: %i" % model.internal_counters[3])
    print("- cost: %f" % cost)
    print("")


def run_calc_cost_and_grad(config, display_cost=True, display_internal_counters=False):
  
    # Create model
    model = create_model(config)
  
    # Load observations
    obs = load_observations(config, model)
  
    # Setup control
    control = load_control(config, model)
    
    # Run gradient test 
    print("=" * 80)
    if display_cost:
        print(" COMPUTE COST AND GRADIENTS")
    else:
        print(" COMPUTE GRADIENTS")
    print("=" * 80)
    
    # Run adjoint model to compute cost and gradient
    model.dtout = -1
    model.disable_stdout = True
    cost, grad = dassflow1d.calc_cost_and_gradients(model, control, obs)
    grad_norm = np.linalg.norm(grad)
    if display_internal_counters:
        print("- Number of overbank flow occurences: %i" % model.internal_counters[0])
        print("  - Left overbank flow occurences : %i" % model.internal_counters[1])
        print("  - Right overbank flow occurences: %i" % model.internal_counters[2])
        print("- Number of overflow (flow above highest elevation) occurences: %i" % model.internal_counters[3])
    if display_cost:
        print("- cost  : %f" % cost)
    print("- grad  : [%f, %f]" % (np.min(grad), np.max(grad)))
    print("- |grad|: %f" % grad_norm)
    print("")


def run_gradient_test(config, iterations=34, delta=1e-3):
  
    # Create model
    model = create_model(config)
  
    # Load observations
    obs = load_observations(config, model)
  
    # Setup control
    control = load_control(config, model)
    
    # Run gradient test 
    print("=" * 80)
    print(" GRADIENT TEST")
    print("=" * 80)
    
    # Run adjoint model to compute cost and gradient
    model.dtout = -1
    model.disable_stdout = True
    cost, grad = dassflow1d.calc_cost_and_gradients(model, control, obs)
    grad_norm = np.linalg.norm(grad)
    print("- cost  : %f" % cost)
    print("- |grad|: %f" % grad_norm)
    print("")
    
    # Loop on epsilon values
    dx = control.x * delta
    dx[np.abs(dx) < 1e-15] = 1.e-8
    alphas = np.zeros(iterations)
    oneminusI_right = np.zeros(iterations)
    oneminusI_left = np.zeros(iterations)
    oneminusI_centered = np.zeros(iterations)
    print(" Alpha        | |1-Iright|   | |1-Icent.|   | |1-Ileft|    |")
    for i in range(iterations):
        
        alpha = 2**(-i)
        alphas[i] = alpha
        
        # Compute "right" cost
        control.x[:] = control.x0 + alpha * dx
        cost_right = dassflow1d.calc_cost(model, control, obs)
        
        # Compute "left" cost
        control.x[:] = control.x0 - alpha * dx
        cost_left = dassflow1d.calc_cost(model, control, obs)
        
        # Compute left, centered and right |1-I|
        oneminusI_right[i] = np.abs(1.0 - (cost_right - cost) / (alpha * np.dot(grad, dx)))
        oneminusI_left[i] = np.abs(1.0 - (cost - cost_left) / (alpha * np.dot(grad, dx)))
        oneminusI_centered[i] = np.abs(1.0 - (cost_right - cost_left) / (2.0 * alpha * np.dot(grad, dx)))
        print(" %12.5e | %12.5e | %12.5e | %12.5e |" % (alpha, oneminusI_right[i], oneminusI_centered[i], oneminusI_left[i]))
        
    plt.plot(alphas, oneminusI_right, 'b-', label="right")
    plt.plot(alphas, oneminusI_left, 'r-', label="left")
    plt.plot(alphas, oneminusI_centered, 'g-', label="centered")
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\left|1-I_\alpha \right|$")
    plt.loglog()
    plt.legend()
    plt.show()
    


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser("Run a Dassflow-1D case")
    parser.add_argument("config_file", type=str, help="Configuration file")
    parser.add_argument("run_type", type=str, default="direct",
                        choices=["direct", "calc_cost", "calc_grad", "calc_cost_and_grad", "gradient_test", "optim"],
                        help="Type of run")
    parser.add_argument("--display-internal-counters", dest="display_internal_counters", action="store_true", 
                        help="Display internal counters values")
    args = parser.parse_args()
    
    # Load configuration
    config = load_configuration(args.config_file)
    
    # Process run
    if args.run_type == "direct":

        # Direct run
        run_direct(config, display_internal_counters=args.display_internal_counters)

    elif args.run_type == "calc_cost":

        # Direct run and cost computation
        run_calc_cost(config, display_internal_counters=args.display_internal_counters)

    elif args.run_type == "calc_grad":

        # Compute gradients and norm
        run_calc_cost_and_grad(config, display_cost=False, display_internal_counters=args.display_internal_counters)

    elif args.run_type == "calc_cost_and_grad":

        # Compute cost, gradients and norm
        run_calc_cost_and_grad(config, display_internal_counters=args.display_internal_counters)

    elif args.run_type == "gradient_test":

        # Gradient test
        run_gradient_test(config)
