import argparse
import geopandas as gpd
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

def auto_filepath(path, candidate_dirs):

    if isinstance(candidate_dirs, str):
        candidate_dirs = [candidate_dirs]
    for candidate_dir in candidate_dirs:
        candidate_file = os.path.join(candidate_dir, path)
        if os.path.isfile(candidate_file):
            return candidate_file
    return path


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

    # # Try with internal format
    # with open(fname, "r") as fp:
    #     content = fp.readlines()

    # # Count number of lines with sharp at the beginning
    # row_index = 0
    # while content[row_index].rstrip()[0:1] == "#":
    #     row_index += 1
    #     if row_index >= len(content):
    #         break
    # print("row_index =", row_index)
     
    if os.path.splitext(fname)[1] == ".csv":
        
        if datestart is None:
            
            data = pd.read_csv(fname, sep=";")
            return data.iloc[:, 0].values, data.iloc[:, 1].values
        
        else:
            
            data = pd.read_csv(fname, sep=",", parse_dates={"datetime": ["date", "time"]})
            t = ((data.iloc[:, 0].dt.tz_localize(None) - np.datetime64(datestart)) / np.timedelta64(1, "s")).values
            
            return t, data.iloc[:, 1].values
            
    else:
        
        raise RuntimeError("Timeseries file must be a CSV file")


def load_provider_timeseries(fname, datestart=None):

    float_keys = ["REFERENCE LONGITUDE", "REFERENCE_LATITUDE"]


    with open(fname, "r") as fp:
        content = fp.readlines()

    # Read metadata
    metadata = {}
    row_index = 0
    while content[row_index].rstrip()[0:1] == "#":
        row_content = content[row_index].rstrip()[1:]
        if "::" in row_content:
            key, value = row_content.split("::")
            if key in float_keys:
                value = float(value)
            metadata[key] = value
        row_index += 1
        if row_index >= len(content):
            break
     
    data = pd.read_csv(fname, sep="\s+", skiprows=row_index, header=None, parse_dates={"datetime": [0,1]})
    t = ((data.loc[:, "datetime"].dt.tz_localize(None) - np.datetime64(datestart)) / np.timedelta64(1, "s")).values
    
    return t, data.loc[:, 2].values, metadata


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
    mesh_fname = auto_filepath(config["model"]["mesh_file"], "input/static_data")
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

    elif config["model"]["strickler_values"]["type"] == "segment_uniform":
        
        # Set segment values
        values = config["model"]["strickler_values"]["values"]
        mesh.set_strickler_fields_segment(values)
        
    elif config["model"]["strickler_values"]["type"] == "fields":
        pass
        #if "strickler_fields" in config["model"]:
            #if config["model"]["strickler_fields"] == "segment":
                #pass
        #else:
            #raise RuntimeError("'strickler_fields' parameter not specified in 'model' section")
    else:
        raise RuntimeError("wrong type of strickler values : %s" % config["model"]["strickler_values"]["type"])

    # print(mesh.strickler_type_code)
    # print(mesh.cs[2].strickler_params)

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
            # print("date_start=", date_start)
            t, y = load_timeseries(bc_config["timeseries_file"], datestart=date_start)
            model.bc[ibc].set_timeseries(t, y)
        elif "provider" in bc_config:
            candidate_file = auto_filepath("BC%04i.csv" % ibc, "input/dynamic_data")
            if os.path.isfile(candidate_file):
                t, y = load_timeseries(candidate_file, datestart=date_start)
                model.bc[ibc].set_timeseries(t, y)
            elif "code_station" in bc_config:
                bc_file = "%s_%s.csv" % (bc_config["id"], bc_config["code_station"])
                candidate_file = auto_filepath(bc_file, "input/dynamic_data")
                if os.path.isfile(candidate_file):
                    t, y = load_timeseries(candidate_file, datestart=date_start)
                    model.bc[ibc].set_timeseries(t, y)
                else:
                    raise RuntimeError("Unable to load values for boundary condition %i" % ibc)
            else:
                raise RuntimeError("Unable to load values for boundary condition %i" % ibc)
    
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
            # TODO ipmlement auto filenames as for BC (see)
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
    if "timestep" in config["model"]:
        model.dt = config["model"]["timestep"]
    elif "dt" in config["model"]:
        model.dt = config["model"]["dt"]
    else:
        model.dt = 300
    
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

    # Count number of stations
    nobs = 0
    for iobs, obs_config in enumerate(config["observations"]):
        if "timeseries_file" in obs_config:
            nobs += 1
        elif "provider" in obs_config:
            if obs_config["provider"] == "hydroweb":
                # TODO
                raise NotImplementedError("Listing of Hydroweb observations files is not implemented yet")
            elif obs_config["provider"] == "schapi":
                if "code_stations" in obs_config:
                    nobs += len(obs_config["code_stations"])
                else:
                    candidates_files = os.path.join("input", "assim_data", "schapi_*.csv")
                    nobs += len(candidates_files)

                # raise NotImplementedError("Listing of SCHAPI observations files is not implemented yet")
            elif obs_config["provider"] == "icesat2":
                # TODO
                raise NotImplementedError("Listing of ICESat-2 observations files is not implemented yet")
            else:
                raise ValueError("Wrong provider: %s" % obs_config["provider"])
        else:
            raise ValueError("Unable to process observation %i. Either provide timeseries file or provider." % iobs)
    print("- Number of stations: %i" % nobs)
    
    # Allocate observations object
    # nobs = len(config["observations"])
    obs = m_obs.Observations(nobs)
    
    # Setup Verdun-sur-Garonne station
    ista = 0
    for iobs, obs_config in enumerate(config["observations"]):
        
        # Load timeseries
        if "timeseries_file" in obs_config:

            tobs, Hobs = load_timeseries(obs_config["timeseries_file"], datestart=date_start)
        
            # Create array of observed H and W, restricted to the simulation window
            Hobs = Hobs[tobs <= model.te]
            tobs = tobs[tobs <= model.te]
            HWobs = np.ones((2, Hobs.size)) * -1e+99
            HWobs[0, :] = Hobs
            
            # Setup station
            obs.stations[iobs].setup(model.msh, tobs, HWobs, indices=obs_config["index"])
            ista += 1

        elif "provider" in obs_config:

            # Load xs.shp file
            if os.path.isfile("input/static_data/xs.shp"):
                xs_metadata = gpd.read_file("input/static_data/xs.shp")
                raise NotImplementedError("Findin cross-section index from xs.shp is not implemented yet.")
            else:
                xs_metadata = None
                # #TODO put an error message
                # xs_index = 2


            if obs_config["provider"] == "hydroweb":
                # TODO
                raise NotImplementedError("Listing of Hydroweb observations files is not implemented yet")
            elif obs_config["provider"] == "schapi":
                if "code_stations" in obs_config:
                    obs_files = [os.path.join("input", "assim_data", "schapi", "schapi_%s.txt" % code) for code in obs_config["code_stations"]]
                else:
                    obs_files = os.path.join("input", "assim_data", "schapi", "schapi_*.txt")
                for obs_file in obs_files:
                    tobs, Hobs, metadata = load_provider_timeseries(obs_file, datestart=date_start)
                    Hobs = Hobs[tobs <= model.te]
                    tobs = tobs[tobs <= model.te]
                    HWobs = np.ones((2, Hobs.size)) * -1e+99
                    HWobs[0, :] = Hobs

                    if len(tobs) == 0:
                        raise RuntimeError("No valid observations found in file: %s" % obs_file)

                    # Setup station
                    if xs_metadata is not None:
                        # TODO Find index using metadata
                        obs.stations[ista].setup(model.msh, tobs, HWobs, indices=xs_index)
                    else:
                        x = metadata["REFERENCE LONGITUDE"]
                        y = metadata["REFERENCE LATITUDE"]
                        obs.stations[ista].setup(model.msh, tobs, HWobs, coords=[x, y])
                    ista += 1

                # # TODO
                # raise NotImplementedError("Listing of SCHAPI observations files is not implemented yet")
            elif obs_config["provider"] == "icesat2":
                # TODO
                raise NotImplementedError("Listing of ICESat-2 observations files is not implemented yet")

            # candidate_file = auto_filepath("%s_MBC%04i.csv" % ibc, "dynamic_data")
            # if os.path.isfile(candidate_file):
            #     t, y = load_timeseries(candidate_file, datestart=date_start)
            #     model.bc[ibc].set_timeseries(t, y)
            # elif "code_station" in bc_config:
            #     bc_file = "%s_%s.csv" % (bc_config["id"], bc_config["code_station"])
            #     candidate_file = auto_filepath(bc_file, "dynamic_data")
            #     if os.path.isfile(candidate_file):
            #         t, y = load_timeseries(candidate_file, datestart=date_start)
            #         model.bc[ibc].set_timeseries(t, y)
            #     else:
            #         raise RuntimeError("Unable to load values for boundary condition %i" % ibc)
            # else:
            #     raise RuntimeError("Unable to load values for boundary condition %i" % ibc)

        
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
    use_becm = False
    becm_sigmas = []
    becm_lengths = []
    for item_config in config["control"]:
        
        if item_config["variable"] == "Q" or item_config["variable"] == "Qin":

            # Add Upstream bc in control
            for ibc in range(0, len(model.bc)):
                if model.bc[ibc].id == "discharge":
                    control.add_bc_in_control(model, ibc)

        elif item_config["variable"] == "K":

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

        # Setup Background Error Covariance Matrix (BECM) parameters
        if "sigma" in item_config:
            use_becm = True
            becm_sigmas.append(item_config["sigma"])
        else:
            becm_sigmas.append(np.nan)
        if "correlation_length" in item_config:
            use_becm = True
            becm_lengths.append(item_config["correlation_length"])
        else:
            becm_lengths.append(np.nan)

    print("Control size: %i" % control.x.size)

    if use_becm is True:
        print("Use Background Error Correlation Matrix")
        control.x0 = control.x[:]
        # x_prior = control.x0.copy()
        B = control.zero_block_matrix()

        for i, item_config in enumerate(config["control"]):

            sigma = becm_sigmas[i]
            length = becm_lengths[i]
            if np.isfinite(sigma):
                if np.isfinite(length):
                    if item_config["variable"] in ["K", "alpha", "beta", "KLOB", "KCH", "KROB", "bathy"]:
                        B.blocks[i].m[:, :] = sigma**2 * control.spatial_correlation_array(mesh, length)
                    elif item_config["variable"][0:3] == "QIN":
                        ibc = int(item_config["variable"][3:])
                        B.blocks[i].m[:, :] = sigma**2 * control.temporal_correlation_array(model.bc[ibc].ts.t, length)
                else:
                    B.blocks[i].m[:, :] = sigma**2 * np.eye(B.blocks[i].m.shape[0])
            elif np.isfinite(length):
                raise ValueError("Cannot use correlation length without sigma")

        control.set_prior_error_covariance_matrixblock(B)
        control.x[:] = 0.0
    print("")
    
    return control


def run_direct(config, display_internal_counters=False):
  
    # Create model
    model = create_model(config)
    print(model)
    print("")
    # print(model.bc[1].ts.t)
    # print(model.bc[1].ts.y)
    # print(model.msh.strickler_type)
    # print(model.msh.cs[2].strickler_params)
    
    # Run direct
    print("=" * 80)
    print(" SIMULATION")
    print("=" * 80)
    
    # Run direct model
    model.run_unsteady()
    # print(model.status)
    # print(model.res.h[:, 0])
    # bathy = model.msh.get_segment_field(0, "bathy", base_cs=False)
    # x = model.msh.get_segment_field(0, "x", base_cs=False)
    # plt.plot(x, bathy, "k-")
    # plt.plot(x, bathy+model.res.h[2:-2, 0], "b-")
    # plt.show()
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

    plt.plot(obs.obs[0, :], "ro")
    plt.plot(obs.est[0, :], "b-")
    plt.show()

    delta = np.mean(obs.obs[0, :]) - np.mean(obs.est[0, :])
    plt.plot(obs.obs[0, :], "ro")
    plt.plot(obs.est[0, :] + delta, "b-")
    plt.show()

    # print(model.res.t)
    # index = np.argmin(5.09400E+05-model.res.t)
    # bathy = model.msh.get_segment_field(0, "bathy", base_cs=False)
    # x = model.msh.get_segment_field(0, "x", base_cs=False)
    # plt.plot(x, bathy, "k-")
    # plt.plot(x, bathy+model.res.h[2:-2, index], "b-")
    # plt.title("t=%.3f" % model.res.t[index])
    # plt.show()



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
    


def run_optim(config, max_iterations=200, feps=1e-3):

    # Create model
    model = create_model(config)
  
    # Load observations
    obs = load_observations(config, model)
  
    # Setup control
    control = load_control(config, model)
    
    # Run gradient test 
    print("=" * 80)
    print(" OPTIMISATION")
    print("=" * 80)
    
    # Run adjoint model to compute cost and gradient
    dtout0 = model.dtout
    model.dtout = -1
    model.disable_stdout = True

    res = scipy.optimize.minimize(dassflow1d.calc_cost_and_gradients_scipy_minimize, control.x, 
                                  args=(model, control, obs), jac=True, tol = 1e-7, 
                                  method='L-BFGS-B', options={"disp": True})
    
    model.dtout = dtout0
    model.run_unsteady()

    valid = np.ravel(np.argwhere(obs.est[0, :] > -99999))
    yobs = obs.obs[0, valid]
    residuals = obs.est[0, valid] - yobs
    rmse = np.sqrt(np.sum(np.ravel(residuals)**2) / obs.obs.shape[1])
    nse = 1.0 - np.sum(np.ravel(residuals)**2) / np.sum((np.ravel(yobs) - np.mean(np.ravel(yobs)))**2)
    plt.plot(obs.obs[0, :], "ro")
    plt.plot(obs.est[0, :], "b-")
    plt.title("valid: %i, rmse=%.2f cm" % (valid.size, rmse))
    plt.show()

    delta = np.mean(obs.obs[0, :]) - np.mean(obs.est[0, :])
    rmse = np.sqrt(np.sum(np.ravel(residuals+delta)**2) / obs.obs.shape[1])
    plt.plot(obs.obs[0, :], "ro")
    plt.plot(obs.est[0, :] + delta, "b-")
    plt.title("valid: %i, rmse=%.2f cm" % (valid.size, rmse))
    plt.show()


    # cost, grad = dassflow1d.calc_cost_and_gradients(model, control, obs)
    # grad_norm = np.linalg.norm(grad)
    # print("- cost  : %f" % cost)
    # print("- |grad|: %f" % grad_norm)
    # print("")
    
    # # Loop on epsilon values
    # dx = control.x * delta
    # dx[np.abs(dx) < 1e-15] = 1.e-8
    # alphas = np.zeros(iterations)
    # oneminusI_right = np.zeros(iterations)
    # oneminusI_left = np.zeros(iterations)
    # oneminusI_centered = np.zeros(iterations)
    # print(" Alpha        | |1-Iright|   | |1-Icent.|   | |1-Ileft|    |")
    # for i in range(iterations):
        
    #     alpha = 2**(-i)
    #     alphas[i] = alpha
        
    #     # Compute "right" cost
    #     control.x[:] = control.x0 + alpha * dx
    #     cost_right = dassflow1d.calc_cost(model, control, obs)
        
    #     # Compute "left" cost
    #     control.x[:] = control.x0 - alpha * dx
    #     cost_left = dassflow1d.calc_cost(model, control, obs)
        
    #     # Compute left, centered and right |1-I|
    #     oneminusI_right[i] = np.abs(1.0 - (cost_right - cost) / (alpha * np.dot(grad, dx)))
    #     oneminusI_left[i] = np.abs(1.0 - (cost - cost_left) / (alpha * np.dot(grad, dx)))
    #     oneminusI_centered[i] = np.abs(1.0 - (cost_right - cost_left) / (2.0 * alpha * np.dot(grad, dx)))
    #     print(" %12.5e | %12.5e | %12.5e | %12.5e |" % (alpha, oneminusI_right[i], oneminusI_centered[i], oneminusI_left[i]))
        
    # plt.plot(alphas, oneminusI_right, 'b-', label="right")
    # plt.plot(alphas, oneminusI_left, 'r-', label="left")
    # plt.plot(alphas, oneminusI_centered, 'g-', label="centered")
    # plt.xlabel(r"$\alpha$")
    # plt.ylabel(r"$\left|1-I_\alpha \right|$")
    # plt.loglog()
    # plt.legend()
    # plt.show()
    print(res)
    

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

    elif args.run_type == "optim":

        # Gradient test
        run_optim(config)
