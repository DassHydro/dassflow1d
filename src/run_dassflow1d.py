import argparse
import geopandas as gpd
import glob
import io
import json
import logging
import os
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
import scipy.optimize
import sys

import dassflow1d
import dassflow1d.m_control as m_control
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_obs as m_obs


__rootdir__ = "/mnt/rundir"


def auto_filepath(path, candidate_dirs):

    logger = logging.getLogger("dassflow1d")

    if isinstance(candidate_dirs, str):
        candidate_dirs = [candidate_dirs]
    for candidate_dir in candidate_dirs:
        candidate_file = os.path.join(candidate_dir, path)
        logger.debug("candidate_file: %s (found:%s)" % (candidate_file, str(os.path.isfile(candidate_file))))
        if os.path.isfile(candidate_file):
            return candidate_file
    return path


def create_logger(debug: bool=False):
    """ Create logger

        Parameters
        ----------
            debug: bool
                True to enable debug logging level
        Return
        ------
            logger: logging.Logger
    """

    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
    if not type(sys.stdout) == io.TextIOWrapper:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger = logging.getLogger("dassflow1d")
    logger.addHandler(handler)

    # Set level
    if debug is True:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    return logger


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
            
            # data = pd.read_csv(fname, sep=",", parse_dates={"datetime": ["date", "time"]})
            data = pd.read_csv(fname, sep=",")
            data["datetime"] = pd.to_datetime(data["date"] + "T" + data["time"])
            t = ((data.loc[:, "datetime"].dt.tz_localize(None) - np.datetime64(datestart)) / np.timedelta64(1, "s")).values
            
            return t, data.iloc[:, 2].values
            
    else:
        
        raise RuntimeError("Timeseries file must be a CSV file")


def load_provider_timeseries(fname, datestart=None):

    str_keys = ["BASIN", "RIVER"]
    int_keys = ["REFERENCE DISTANCE (km)", "ID"]
    float_keys = ["REFERENCE LONGITUDE", "REFERENCE LATITUDE"]

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
            elif key in str_keys:
                value = value.strip()
            elif key in int_keys:
                value = int(value)
            else:
                value = value.strip()
            metadata[key] = value
        row_index += 1
        if row_index >= len(content):
            break
     
    data = pd.read_csv(fname, sep="\s+", skiprows=row_index, header=None)
    data["datetime"] = pd.to_datetime(data[0] + "T" + data[1])
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

    logger = logging.getLogger("dassflow1d")

    logger.info("=" * 80)
    logger.info(" SETUP MODEL")
    logger.info("=" * 80)

    static_data_dir = os.path.join(__rootdir__, "input", "static_data")
    dynamic_data_dir = os.path.join(__rootdir__, "input", "dynamic_data")
    output_dir = os.path.join(__rootdir__, "output")
    
    # Load mesh
    mesh_fname = auto_filepath(config["model"]["mesh_file"], static_data_dir)
    logger.debug("mesh_fname=", mesh_fname)
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
            candidate_file = auto_filepath("BC%04i.csv" % ibc, dynamic_data_dir)
            if os.path.isfile(candidate_file):
                t, y = load_timeseries(candidate_file, datestart=date_start)
                model.bc[ibc].set_timeseries(t, y)
            elif "code_station" in bc_config:
                bc_file = "%s_%s.csv" % (bc_config["id"], bc_config["code_station"])
                candidate_file = auto_filepath(bc_file, dynamic_data_dir)
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

    if "output_format" in config["model"]:
        output_format = config["model"]["output_format"]
        if output_format not in ["csv", "netcdf"]:
            raise ValueError("Output format must be either 'csv' or 'netcdf'")
    else:
        output_format = "netcdf"
    model.output_format = output_format
    
    # Set output file
    if "output_file" in config["model"]:
        model.output_file = config["model"]["output_file"]
    else:
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        if output_format == "csv":
            model.output_file = os.path.join(output_dir, "results.csv")
        else:
            model.output_file = os.path.join(output_dir, "dassflow1d_out.nc")

    if output_format == "netcdf":
        model.output_nc_file = model.output_file.decode("utf-8")
        model.output_file = ""
    
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

    logger = logging.getLogger("dassflow1d")

    logger.info("=" * 80)
    logger.info(" SETUP OBSERVATIONS STATIONS")
    logger.info("=" * 80)

    static_data_dir = os.path.join(__rootdir__, "input", "static_data")
    assim_data_dir = os.path.join(__rootdir__, "input", "assim_data")
    # dynamic_data_dir = os.path.join(__rootdir__, "input", "dynamic_data")
    # output_dir = os.path.join(__rootdir__, "output")
    
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

                logger.debug("Search for hydroweb files")
                candidates_files = glob.glob(os.path.join(assim_data_dir, "hydroweb", "hydroweb_*.txt"))
                logger.debug("- Number of hydroweb files found: %i" % len(candidates_files))
                nobs += len(candidates_files)

            elif obs_config["provider"] == "schapi":

                if "code_stations" in obs_config:
                    logger.debug("Add %i schapi stations from codes defined in the config file")
                    nobs += len(obs_config["code_stations"])
                else:
                    logger.debug("Search for schapi files")
                    candidates_files = os.path.join(assim_data_dir, "schapi_*.csv")
                    logger.debug("- Number of schapi files found: %i" % len(candidates_files))
                    nobs += len(candidates_files)

                # raise NotImplementedError("Listing of SCHAPI observations files is not implemented yet")
            elif obs_config["provider"] == "icesat2":
                # TODO
                raise NotImplementedError("Listing of ICESat-2 observations files is not implemented yet")
            else:
                raise ValueError("Wrong provider: %s" % obs_config["provider"])
        else:
            raise ValueError("Unable to process observation %i. Either provide timeseries file or provider." % iobs)
    logger.info("- Number of stations: %i" % nobs)
    
    # Allocate observations object
    # nobs = len(config["observations"])
    obs = m_obs.Observations(nobs)
    
    # Setup observations station
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
            if os.path.isfile(os.path.join(static_data_dir, "xs.shp")):
                xs_metadata = gpd.read_file(os.path.join(static_data_dir, "xs.shp"))
                raise NotImplementedError("Findin cross-section index from xs.shp is not implemented yet.")
            elif os.path.isfile(os.path.join(static_data_dir, "xs_nodes.shp")):
                logger.debug("Use xs<->nodes assocation file for provider %s" % obs_config["provider"])
                xs_metadata = gpd.read_file(os.path.join(static_data_dir, "xs_nodes.shp"))
            else:
                xs_metadata = None

            if obs_config["provider"] == "hydroweb":

                logger.debug("Load hydroweb observations files")
                obs_files = glob.glob(os.path.join(assim_data_dir, "hydroweb", "hydroweb_*.txt"))

                for obs_file in obs_files:

                    logger.debug("- Load hydroweb observations file: %s" % os.path.basename(obs_file))
                    tobs, Hobs, metadata = load_provider_timeseries(obs_file, datestart=date_start)
                    if np.any(tobs <= model.te):
                        Hobs = Hobs[tobs <= model.te]
                        tobs = tobs[tobs <= model.te]
                    HWobs = np.ones((2, Hobs.size)) * -1e+99
                    HWobs[0, :] = Hobs

                    if len(tobs) == 0:
                        raise RuntimeError("No valid observations found in file: %s" % obs_file)
                    logger.debug("  - Number of observations: %i" % len(tobs))

                    # Setup station
                    xs_index = None
                    if xs_metadata is not None:

                        node_id = metadata["ID"]
                        candidates = xs_metadata[xs_metadata["node_id"] == node_id]
                        if candidates.index.size == 1:
                            xs_index = candidates["xs_idx"].values[0]

                    if xs_index is not None:
                        logger.debug("  - Station index: %i (ncs=%i)" % (xs_index, model.msh.ncs))
                        obs.stations[ista].setup(model.msh, tobs, HWobs, indices=xs_index)
                    else:
                        x = metadata["REFERENCE LONGITUDE"]
                        y = metadata["REFERENCE LATITUDE"]
                        logger.debug("  - Station coordinates: %f,%f" % (x, y))
                        obs.stations[ista].setup(model.msh, tobs, HWobs, coords=[x, y])
                        logger.debug("  - Station index: %i" % obs.stations[ista].ics[0])
                    ista += 1
            
            elif obs_config["provider"] == "schapi":

                logger.debug("Load schapi observations files")
                if "code_stations" in obs_config:
                    obs_files = [os.path.join(assim_data_dir, "schapi", "schapi_%s.txt" % code) for code in obs_config["code_stations"]]
                else:
                    obs_files = glob.glob(os.path.join(assim_data_dir, "schapi", "schapi_*.txt"))
                for obs_file in obs_files:
                    logger.debug("- Load schapi observations file: %s" % os.path.basename(obs_file))
                    tobs, Hobs, metadata = load_provider_timeseries(obs_file, datestart=date_start)
                    Hobs = Hobs[tobs <= model.te]
                    tobs = tobs[tobs <= model.te]
                    HWobs = np.ones((2, Hobs.size)) * -1e+99
                    HWobs[0, :] = Hobs

                    if len(tobs) == 0:
                        raise RuntimeError("No valid observations found in file: %s" % obs_file)
                    logger.debug("  - Number of observations: %i" % len(tobs))

                    # Setup station
                    # TODO write a common method (is identical to hydroweb association)
                    xs_index = None
                    if xs_metadata is not None:

                        node_id = metadata["ID"]
                        candidates = xs_metadata[xs_metadata["node_id"] == node_id]
                        if candidates.index.size == 1:
                            xs_index = candidates["xs_idx"].values[0]
                        # else:
                        #     xs_index = None
                        #     # raise RuntimeError("No XS associated to node %i" % node_id)
                        # # elif candidates.index.size > 1:
                        # #     raise RuntimeError("Multiple XS associated to node %i" % node_id)
                        # # xs_index = candidates["xs_idx"].values[0]
                        # obs.stations[ista].setup(model.msh, tobs, HWobs, indices=xs_index)

                    if xs_index is not None:
                        logger.debug("  - Station index: %i (ncs=%i)" % (xs_index, model.msh.ncs))
                        obs.stations[ista].setup(model.msh, tobs, HWobs, indices=xs_index)
                    else:
                        x = metadata["REFERENCE LONGITUDE"]
                        y = metadata["REFERENCE LATITUDE"]
                        logger.debug("  - Station coordinates: %f,%f" % (x, y))
                        obs.stations[ista].setup(model.msh, tobs, HWobs, coords=[x, y])
                        logger.debug("  - Station index: %i" % obs.stations[ista].ics[0])
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

        
    logger.info("")
    
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

    logger = logging.getLogger("dassflow1d")

    logger.info("=" * 80)
    logger.info(" SETUP CONTROL")
    logger.info("=" * 80)

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

    logger.info("Control size: %i" % control.x.size)

    if use_becm is True:
        logger.info("Use Background Error Correlation Matrix")
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
    logger.info("")
    
    return control

def write_netcdf_results(model, config):

    logger = logging.getLogger("dassflow1d")
    logging.info("Export results to NetCDF file")

    # Compute dates array
    date_start = np.datetime64(config["model"]["date_start"])
    nt = len(model.res.t)
    dates = np.array([date_start + np.timedelta64(int(x * model.dtout), "s") for x in range(0, nt)])

    # Retrieve coordinates of cross-sections
    mesh = model.msh
    coords = np.zeros(shape=(0,2))
    for iseg in range(0, mesh.nseg):
        coords = np.concatenate((coords, mesh.get_segment_field(iseg, "coords")), axis=0)
    nx = coords.shape[0]

    # Compute indices of real cross-sections in res array
    indices = []
    for iseg in range(0, mesh.nseg):
        for ics in range(mesh.seg[iseg].first_cs-1, mesh.seg[iseg].last_cs):
            if mesh.cs[ics].ibase > 0:
                # print(ics, mesh.cs[ics].ibase)
                indices.append(ics-1)
        #         choice = input()
        # coords = np.concatenate((coords, mesh.get_segment_field(iseg, "coord")), axis=0)

    # Create netCDF output file
    logger.info("- Create NetCDF output file: %s" % model.output_nc_file)
    dataset = nc.Dataset(model.output_nc_file, "w", format="NETCDF4")
    dataset.createDimension("nx", nx)
    dataset.createDimension("nt", nt)
    time = dataset.createVariable("time", "f4", ("nt",))
    time.setncattr("units", "seconds since 2000-01-01 00:00:00")
    time[:] = (dates - np.datetime64("2000-01-01 00:00:00")) / np.timedelta64(1, "s")
    x = dataset.createVariable("x", "f4", ("nx",))
    x[:] = coords[:, 0]
    y = dataset.createVariable("y", "f4", ("nx",))
    y[:] = coords[:, 1]
    h = dataset.createVariable("h", "f4", ("nt", "nx"))
    h.setncattr("units", "m3/s")
    offset = 0
    h[:, :] = model.res.h[indices, :].T
    Q = dataset.createVariable("Q", "f4", ("nt", "nx"))
    Q.setncattr("units", "m3/s")
    Q[:, :] = model.res.q[indices, :].T
    dataset.close()





def run_direct(config, display_internal_counters=False):

    logger = logging.getLogger("dassflow1d")
  
    # Create model
    model = create_model(config)
    # print(model.bc[1].ts.t)
    # print(model.bc[1].ts.y)
    # print(model.msh.strickler_type)
    # print(model.msh.cs[2].strickler_params)
    
    # Run direct
    logger.info("=" * 80)
    logger.info(" SIMULATION")
    logger.info("=" * 80)
    
    # Run direct model
    model.run_unsteady()

    if model.output_format == "netcdf":
        write_netcdf_results(model, config)
    # print(model.status)
    # print(model.res.h[:, 0])
    # bathy = model.msh.get_segment_field(0, "bathy", base_cs=False)
    # x = model.msh.get_segment_field(0, "x", base_cs=False)
    # plt.plot(x, bathy, "k-")
    # plt.plot(x, bathy+model.res.h[2:-2, 0], "b-")
    # plt.show()
    if display_internal_counters:
        logger.debug("- Number of overbank flow occurences: %i" % model.internal_counters[0])
        logger.debug("  - Left overbank flow occurences : %i" % model.internal_counters[1])
        logger.debug("  - Right overbank flow occurences: %i" % model.internal_counters[2])
        logger.debug("- Number of overflow (flow above highest elevation) occurences: %i" % model.internal_counters[3])
    logger.info("")


def run_calc_cost(config, display_internal_counters=False):

    logger = logging.getLogger("dassflow1d")
  
    # Create model
    model = create_model(config)
  
    # Load observations
    obs = load_observations(config, model)
  
    # Setup dummy control
    control = m_control.Control()
    
    # Run direct 
    logger.info("=" * 80)
    logger.info(" SIMULATION")
    logger.info("=" * 80)
    
    # Run direct model to compute cost
    cost = dassflow1d.calc_cost(model, control, obs)
    if display_internal_counters:
        logger.debug("- Number of overbank flow occurences: %i" % model.internal_counters[0])
        logger.debug("  - Left overbank flow occurences : %i" % model.internal_counters[1])
        logger.debug("  - Right overbank flow occurences: %i" % model.internal_counters[2])
        logger.debug("- Number of overflow (flow above highest elevation) occurences: %i" % model.internal_counters[3])
    logger.info("- cost: %f" % cost)
    logger.info("")

    # plt.plot(obs.obs[0, :], "ro")
    # plt.plot(obs.est[0, :], "b-")
    # plt.show()

    # delta = np.mean(obs.obs[0, :]) - np.mean(obs.est[0, :])
    # plt.plot(obs.obs[0, :], "ro")
    # plt.plot(obs.est[0, :] + delta, "b-")
    # plt.show()

    # print(model.res.t)
    # index = np.argmin(5.09400E+05-model.res.t)
    # bathy = model.msh.get_segment_field(0, "bathy", base_cs=False)
    # x = model.msh.get_segment_field(0, "x", base_cs=False)
    # plt.plot(x, bathy, "k-")
    # plt.plot(x, bathy+model.res.h[2:-2, index], "b-")
    # plt.title("t=%.3f" % model.res.t[index])
    # plt.show()



def run_calc_cost_and_grad(config, display_cost=True, display_internal_counters=False):

    logger = logging.getLogger("dassflow1d")

    # Create model
    model = create_model(config)
  
    # Load observations
    obs = load_observations(config, model)
  
    # Setup control
    control = load_control(config, model)
    
    # Run gradient test 
    logger.info("=" * 80)
    if display_cost:
        logger.info(" COMPUTE COST AND GRADIENTS")
    else:
        logger.info(" COMPUTE GRADIENTS")
    logger.info("=" * 80)
    
    # Run adjoint model to compute cost and gradient
    model.dtout = -1
    model.disable_stdout = True
    cost, grad = dassflow1d.calc_cost_and_gradients(model, control, obs)
    grad_norm = np.linalg.norm(grad)
    if display_internal_counters:
        logger.debug("- Number of overbank flow occurences: %i" % model.internal_counters[0])
        logger.debug("  - Left overbank flow occurences : %i" % model.internal_counters[1])
        logger.debug("  - Right overbank flow occurences: %i" % model.internal_counters[2])
        logger.debug("- Number of overflow (flow above highest elevation) occurences: %i" % model.internal_counters[3])
    if display_cost:
        logger.info("- cost  : %f" % cost)
    logger.info("- grad  : [%f, %f]" % (np.min(grad), np.max(grad)))
    logger.info("- |grad|: %f" % grad_norm)
    logger.info("")


def run_gradient_test(config, iterations=34, delta=1e-3):

    logger = logging.getLogger("dassflow1d")

    # Create model
    model = create_model(config)
  
    # Load observations
    obs = load_observations(config, model)
  
    # Setup control
    control = load_control(config, model)
    
    # Run gradient test 
    logger.info("=" * 80)
    logger.info(" GRADIENT TEST")
    logger.info("=" * 80)
    
    # Run adjoint model to compute cost and gradient
    model.dtout = -1
    model.disable_stdout = True
    cost, grad = dassflow1d.calc_cost_and_gradients(model, control, obs)
    grad_norm = np.linalg.norm(grad)
    logger.info("- cost  : %f" % cost)
    logger.info("- |grad|: %f" % grad_norm)
    logger.info("")
    
    # Loop on epsilon values
    dx = control.x * delta
    dx[np.abs(dx) < 1e-15] = 1.e-8
    alphas = np.zeros(iterations)
    oneminusI_right = np.zeros(iterations)
    oneminusI_left = np.zeros(iterations)
    oneminusI_centered = np.zeros(iterations)
    logger.info(" Alpha        | |1-Iright|   | |1-Icent.|   | |1-Ileft|    |")
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
        logger.info(" %12.5e | %12.5e | %12.5e | %12.5e |" % (alpha, oneminusI_right[i], oneminusI_centered[i], oneminusI_left[i]))
        
    # plt.plot(alphas, oneminusI_right, 'b-', label="right")
    # plt.plot(alphas, oneminusI_left, 'r-', label="left")
    # plt.plot(alphas, oneminusI_centered, 'g-', label="centered")
    # plt.xlabel(r"$\alpha$")
    # plt.ylabel(r"$\left|1-I_\alpha \right|$")
    # plt.loglog()
    # plt.legend()
    # plt.show()
    


def run_optim(config, max_iterations=200, feps=1e-3):

    logger = logging.getLogger("dassflow1d")

    # Create model
    model = create_model(config)
  
    # Load observations
    obs = load_observations(config, model)
  
    # Setup control
    control = load_control(config, model)
    
    # Run gradient test 
    logger.info("=" * 80)
    logger.info(" OPTIMISATION")
    logger.info("=" * 80)
    
    # Run adjoint model to compute cost and gradient
    dtout0 = model.dtout
    model.dtout = -1
    model.disable_stdout = True

    def callback_minimize(x):
        logger.info("Iteration, cost= %12.5e, [proj g]= %12.5e" % (model.__last_cost__, np.linalg.norm(model.__last_grad__)))

    res = scipy.optimize.minimize(dassflow1d.calc_cost_and_gradients_scipy_minimize, control.x, 
                                  args=(model, control, obs), jac=True, tol = 1e-7, 
                                  method='L-BFGS-B', options={"disp": False},
                                  callback=callback_minimize)
    
    if res.success is True:
        logger.info("Optimisation is successful")
    else:
        logger.info("Optimisation has failed: %s" % res.message)
    logger.info("- Number of iterations: %s" % res.nit)
    logger.info("- Number of cost evaluation: %s" % res.nfev)
    
    model.dtout = dtout0
    model.print_progress = False
    logger.info("Compute results with optimal control")
    model.run_unsteady()

    if model.output_format == "netcdf":
        write_netcdf_results(model, config)


    # valid = np.ravel(np.argwhere(obs.est[0, :] > -99999))
    # yobs = obs.obs[0, valid]
    # residuals = obs.est[0, valid] - yobs
    # rmse = np.sqrt(np.sum(np.ravel(residuals)**2) / obs.obs.shape[1])
    # nse = 1.0 - np.sum(np.ravel(residuals)**2) / np.sum((np.ravel(yobs) - np.mean(np.ravel(yobs)))**2)
    # plt.plot(obs.obs[0, :], "ro")
    # plt.plot(obs.est[0, :], "b-")
    # plt.title("valid: %i, rmse=%.2f cm" % (valid.size, rmse))
    # plt.show()

    # delta = np.mean(obs.obs[0, :]) - np.mean(obs.est[0, :])
    # rmse = np.sqrt(np.sum(np.ravel(residuals+delta)**2) / obs.obs.shape[1])
    # plt.plot(obs.obs[0, :], "ro")
    # plt.plot(obs.est[0, :] + delta, "b-")
    # plt.title("valid: %i, rmse=%.2f cm" % (valid.size, rmse))
    # plt.show()


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
    # print(res)
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser("Run a Dassflow-1D case")
    parser.add_argument("config_file", type=str, help="Configuration file")
    parser.add_argument("run_type", type=str, default="direct",
                        choices=["direct", "calc_cost", "calc_grad", "calc_cost_and_grad", "gradient_test", "optim"],
                        help="Type of run")
    parser.add_argument("--display-internal-counters", dest="display_internal_counters", action="store_true", 
                        help="Display internal counters values")
    parser.add_argument("--debug", action="store_true", 
                        help="Enable debug log")
    args = parser.parse_args()

    # Create Logger
    logger = create_logger(args.debug)
    
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

        # Optim/Assimilation
        run_optim(config)
