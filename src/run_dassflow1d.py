import argparse
import geopandas as gpd
import glob
import io
import json
import logging
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import netCDF4 as nc
import numpy as np
import pandas as pd
import scipy.optimize
import sys
import warnings
# warnings.filterwarnings(action="error")

# import matplotlib as mpl
# mpl.rcParams['text.usetex'] = False

import dassflow1d
import dassflow1d.m_control as m_control
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_obs as m_obs


from datetime import datetime, timedelta
from dassflow1d.post.results import plot_BZQW



__rootdir__ = "/mnt/rundir"

def str2bool(value):

    if isinstance(value, bool):
        return value
    elif value == "true" or value == "True":
        return True
    elif value == "false" or value == "False":
        return False
    else:
        raise ValueError(f"Value '{value}' cannot converted to boolean")
    
def validated_date_string(value):

    try:
        date = np.datetime64(value)
    except:
        raise ValueError(f"Value '{value}' cannot be converted to datetime64")
    return date

def get_from_config(config, section, key, value_type=str, required=True, default=None, deprecated_keys=[]):
    
    # Get config (sub)section
    sections = section.split("/")
    base_section_name = "[config-root]"
    base_section = config
    for subsection in sections:
        if subsection in base_section:
            base_section = base_section[subsection]
            base_section_name += "/%s" % subsection
        elif subsection.replace("-", "_") in base_section:
            logger = logging.getLogger("dassflow1d")
            logger.warning("Subsection '%s' is deprecated in section %s. Replace with '%s'" % (subsection.replace("_", "-"), base_section_name, subsection))
            base_section = base_section[subsection.replace("-", "_")]
            base_section_name += "/%s" % subsection
        elif required:
            raise RuntimeError("subsection %s not found in section %s" % (subsection, base_section_name))
        else:
            return default

    if key in base_section:
        return value_type(base_section[key])
    else:
        for deprecated_key in deprecated_keys:
            if deprecated_key in base_section:
                logger = logging.getLogger("dassflow1d")
                logger.warning("'%s' is deprecated in section %s. Replace with '%s'" % (deprecated_key, base_section_name, key))
                return value_type(base_section[deprecated_key])
        if required:
            raise RuntimeError("key %s not found in section %s" % (key, base_section_name))
        else:
            return default

def get_from_config(config, section, key, value_type=str, required=True, default=None, deprecated_keys=[]):
    
    # Get config (sub)section
    sections = section.split("/")
    base_section_name = "[config-root]"
    base_section = config
    for subsection in sections:
        if subsection in base_section:
            base_section = base_section[subsection]
            base_section_name += "/%s" % subsection
        elif subsection.replace("-", "_") in base_section:
            logger = logging.getLogger("dassflow1d")
            logger.warning("Subsection '%s' is deprecated in section %s. Replace with '%s'" % (subsection.replace("_", "-"), base_section_name, subsection))
            base_section = base_section[subsection.replace("-", "_")]
            base_section_name += "/%s" % subsection
        elif required:
            raise RuntimeError("subsection %s not found in section %s" % (subsection, base_section_name))
        else:
            return default

    if value_type == "index" and isinstance(key, int):
        if required and len(base_section) <= key:
            raise RuntimeError("index %i not found in section %s" % (key, base_section_name))
        return base_section[key]
    if key in base_section:
        return value_type(base_section[key])
    else:
        for deprecated_key in deprecated_keys:
            if deprecated_key in base_section:
                logger = logging.getLogger("dassflow1d")
                logger.warning("'%s' is deprecated in section %s. Replace with '%s'" % (deprecated_key, base_section_name, key))
                return value_type(base_section[deprecated_key])
        if required:
            raise RuntimeError("key %s not found in section %s" % (key, base_section_name))
        else:
            return default


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


def load_timeseries(fname, datestart=None, sep=None):

    if os.path.splitext(fname)[1] == ".csv":
        
        if datestart is None:
            
            data = pd.read_csv(fname, sep=";")
            return data.iloc[:, 0].values, data.iloc[:, 1].values
        
        else:
            
            # data = pd.read_csv(fname, sep=",", parse_dates={"datetime": ["date", "time"]})
            if sep is None:
                sep = ","
            data = pd.read_csv(fname, sep=sep)
            # print(data)
            varcolumn = data.columns.size-1
            if data.columns[0] == "date" and data.columns[1] =="time":
                data["datetime"] = pd.to_datetime(data["date"] + "T" + data["time"])
                varcolumn = 2
            elif data.columns[0] == "t":
                data["datetime"] = [np.datetime64(datestart) + np.timedelta64(int(dt), "s") for dt in data["t"].values]
                varcolumn = 1
            else:
                data["datetime"] = pd.to_datetime(data["date"])
                varcolumn = 1
            t = ((data.loc[:, "datetime"].dt.tz_localize(None) - np.datetime64(datestart)) / np.timedelta64(1, "s")).values
            
            return t, data.iloc[:, varcolumn].values
            
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
    logger.info("Load mesh")
    mesh_fname = get_from_config(config, section="model", 
                                 key="mesh-file", required=True, 
                                 value_type=str, default=None,
                                 deprecated_keys=["mesh_file"])
    mesh_fname = auto_filepath(mesh_fname, static_data_dir)
    logger.debug("Mesh file: %s" % mesh_fname)
    mesh = dassflow1d.read_mesh(mesh_fname)
    logger.info("Mesh has %i cross-sections in %i segments" % (mesh.ncs, mesh.nseg))
    
    # Resample mesh if requested
    mesh_resample = get_from_config(config, section="model", 
                                    key="mesh-resample", required=False,
                                    value_type=float, default=None,
                                    deprecated_keys=["mesh_resample"])
    if mesh_resample is not None:
        mesh.resample(mesh_resample)
        logger.info("Resample mesh at dx: %.3f m" % mesh_resample)
        logger.info("Resampled mesh has %i cross-sections" % mesh.ncs)
        
    # Configure Strickler type
    strickler_type = get_from_config(config, section="model", 
                                     key="strickler-type", required=False, 
                                     value_type=str, default="constant", 
                                     deprecated_keys=["strickler_type"])
    mesh.set_strickler_type(strickler_type)
        
        
    # Setup strickler values
    strickler_values_type = get_from_config(config, section="model/strickler-values", 
                                            key="type", required=True, 
                                            value_type=str, default="constant", 
                                            deprecated_keys=["strickler_type"])
    if strickler_values_type == "uniform":
        
        # Set uniform values
        values = get_from_config(config, section="model/strickler-values", 
                                 key="values", required=True, 
                                 value_type=np.array, default=None)
        # values = config["model"]["strickler_values"]["values"]
        mesh.set_uniform_strickler_parameters(values)

    elif strickler_values_type.replace("_", "-") == "segment-uniform":

        if strickler_values_type == "segment_uniform":
            logger.warning("Strikler values type 'segment_uniform' is deprecated. Replace with 'segment-uniform'")
        
        # Set segment values
        values = get_from_config(config, section="model/strickler-values", 
                                 key="values", required=True, 
                                 value_type=np.array, default=None)
        # values = np.array(config["model"]["strickler_values"]["values"])
        if values.ndim == 1:
            values = np.repeat(values.reshape((values.size, 1)), mesh.nseg, axis=1)
        mesh.set_strickler_fields_segment(values)
        
    elif strickler_values_type == "fields":
        pass
        #if "strickler_fields" in config["model"]:
            #if config["model"]["strickler_fields"] == "segment":
                #pass
        #else:
            #raise RuntimeError("'strickler_fields' parameter not specified in 'model' section")
    else:
        raise RuntimeError("wrong type of strickler values : %s" % strickler_values_type)

    # print(mesh.strickler_type_code)
    # print(mesh.cs[2].strickler_params)

    # Create model from mesh
    logger.info("Init model")
    model = m_sw_mono.Model(mesh)
    
    # Retrieve simulation dates or duration (in seconds)
    date_start = get_from_config(config, section="model", 
                                 key="start-date", required=True, 
                                 value_type=validated_date_string, 
                                 deprecated_keys=["date_start"])
    date_end = get_from_config(config, section="model", 
                               key="end-date", required=False, 
                               value_type=validated_date_string, default=None,
                               deprecated_keys=["date_end"])
    duration = get_from_config(config, section="model", 
                               key="duration", required=False, 
                               value_type=float, default=None)
    # if "date-start" in config["model"]:
    #     date_start = config["model"]["date-start"]
    # elif "date_start" in config["model"]:
    #     # DEPRECATION WARNING
    #     logger.warning("'date_start in model section is deprecated. It must be replaced by 'date-start'")
    #     date_start = config["model"]["date_start"]
    # else:
    #     date_start = None
    # if "date_end" in config["model"]:
    #     date_end = config["model"]["date_end"]
    # else:
    #     date_end = None
    # if "duration" in config["model"]:
    #     duration = config["model"]["duration"]
    # else:
    #     duration = None

    # Set start time
    model.ts = 0
    
    # Set end time
    if duration is None:
        model.te = (np.datetime64(date_end) - np.datetime64(date_start)) / np.timedelta64(1, "s")
    else:
        model.te = duration
    
    # Set minimal depth
    model.heps = get_from_config(config, section="model", 
                                 key="heps", required=False, 
                                 value_type=float, default=0.01)
    # if "heps" in config["model"]:
    #     model.heps = config["model"]["heps"]
    # else:
    #     model.heps = 0.01
    
    # Set minimal discharge
    model.heps = get_from_config(config, section="model", 
                                 key="qeps", required=False, 
                                 value_type=float, default=0.1)
    # if "qeps" in config["model"]:
    #     model.qeps = config["model"]["qeps"]
    # else:
    #     model.qeps = 0.1

    # Load coupling mode
    coupling_mode = get_from_config(config, section="model/coupling", 
                                    key="mode", required=False, 
                                    value_type=str, default="bc+qlat")

    # Load coupling matrix if provided
    coupling_table = get_from_config(config, section="model/coupling", 
                                     key="table-file", required=False, 
                                     value_type=str, default=None,
                                     deprecated_keys=["matrix-file", "matrix_file"])
    if coupling_table is not None:
        table_file = auto_filepath(coupling_table, static_data_dir)
        coupling_table = pd.read_csv(table_file)
        coupling_bc = coupling_table[coupling_table["ics"] == 0].sort_values(by="segment")
    # coupling_table = None
    # if "coupling" in config["model"]:
    #     if "matrix-file" in config["model"]["coupling"]:
    #         matrix_file = auto_filepath(config["model"]["coupling"]["matrix-file"], static_data_dir)
    #         coupling_table = pd.read_csv(matrix_file)
    #         coupling_bc = coupling_table[coupling_table["ics"] == 0].sort_values(by="segment")
    
    # Setup boundary conditions
    ibc_config = 0
    ibc_coupling = 0
    # for ibc, bc_config in enumerate(config["model"]["boundary_conditions"]):
    logger.info("Setup boundary conditions")
    for ibc in range(0, len(model.bc)):
        if coupling_table is not None:
            
            if ibc_coupling < coupling_bc.index.size:

                logger.info("BC %i from %i-th coupling BC" % (ibc, ibc_coupling))

                model.bc[ibc].id = "discharge"
                data_file = auto_filepath("BC%04i.csv" % (ibc_coupling+1), dynamic_data_dir)
                if not os.path.isfile(data_file):
                    logger.error("Unable to find coupling BC file: %s" % data_file)
                    raise IOError("Unable to find coupling BC file: %s" % data_file)
                t, y = load_timeseries(data_file, datestart=date_start)
                y = y[np.logical_and(t >= model.ts, t <= model.te)]
                t = t[np.logical_and(t >= model.ts, t <= model.te)]
                # plt.plot(t, y)
                # plt.title("BC%04i [COUPLING]" % ibc)
                # plt.savefig("plot/direct/debug/BC%03i.png" % ibc)
                # plt.close(plt.gcf())
                model.bc[ibc].set_timeseries(t, np.maximum(model.qeps, y))
                ibc_coupling += 1
                continue

        logger.info("BC %i from config 'boundary-conditions'" % ibc)
        # bc_config = config["model"]["boundary_conditions"][ibc_config]
        bc_config = get_from_config(config, section="model/boundary-conditions", 
                                    key=ibc_config, required=True, 
                                    value_type="index")
        model.bc[ibc].id = bc_config["id"] 
        logger.info("- Type: %s" % bc_config["id"])
        if bc_config["id"] == "normal_depth":
            logger.info("  - Downstream slope: %.3f m/km" % (model.msh.cs[model.msh.seg[model.msh.nseg-1].last_cs-1].slope * 1000.0))
        if "timeseries_file" in bc_config:
            # print("date_start=", date_start)
            # t, y = load_timeseries(bc_config["timeseries_file"], datestart=date_start)
            if "separator" in bc_config:
                sep = bc_config["separator"]
            else:
                sep = None
            t, y = load_timeseries(bc_config["timeseries_file"], datestart=date_start, sep=sep)
            if "correction" in bc_config:
                correction = bc_config["correction"]
                y += correction
            # print("timeseries_file:t=", t)
            # print("timeseries_file:y=", y)
            # choice = input()
            # plt.plot(t, y)
            # plt.title("BC%04i [TS FILE]" % ibc)
            # plt.savefig("plot/direct/debug/BC%03i.png" % ibc)
            # plt.close(plt.gcf())
            model.bc[ibc].set_timeseries(t, y)
        elif "provider" in bc_config:
            candidate_file = auto_filepath("BC%04i.csv" % ibc, dynamic_data_dir)
            if os.path.isfile(candidate_file):
                t, y = load_timeseries(candidate_file, datestart=date_start)
                # plt.plot(t, y)
                # plt.title("BC%04i [PROVIDER-1]" % ibc)
                # plt.savefig("plot/direct/debug/BC%03i.png" % ibc)
                # plt.close(plt.gcf())
                model.bc[ibc].set_timeseries(t, y)
            elif "code_station" in bc_config:
                bc_file = "%s_%s.csv" % (bc_config["id"], bc_config["code_station"])
                candidate_file = auto_filepath(bc_file, dynamic_data_dir)
                if os.path.isfile(candidate_file):
                    t, y = load_timeseries(candidate_file, datestart=date_start)
                    # plt.plot(t, y)
                    # plt.title("BC%04i [PROVIDER-2]" % ibc)
                    # plt.savefig("plot/direct/debug/BC%03i.png" % ibc)
                    # plt.close(plt.gcf())
                    model.bc[ibc].set_timeseries(t, y)
                else:
                    raise RuntimeError("Unable to load values for boundary condition %i" % ibc)
            else:
                raise RuntimeError("Unable to load values for boundary condition %i" % ibc)
            
        ibc_config += 1
            
    # for ibc in range(0, len(model.bc)):
    #     print("BC:", model.bc[ibc].id)
    # choice = input()
    
    # Setup inflow conditions
    if coupling_table is not None and coupling_mode in ["bc+qlat", "qlat"]:

        coupling_inflow = coupling_table[coupling_table["ics"] > 0]

        # Load coupling data
        if "data-file" in config["model"]["coupling"]:
            data_file = config["model"]["coupling"]["data-file"]
        else:
            data_file = "lateral_inflows.nc"
        data_file = auto_filepath(data_file, dynamic_data_dir)
        if not os.path.isfile(data_file):
            logger.error("Unable to find coupling data file: %s" % data_file)
            raise IOError("Unable to find coupling data file: %s" % data_file)
        data = nc.Dataset(data_file, "r")
        time = data.variables["time"][:]
        dates = np.array([np.datetime64("2000-01-01") + np.timedelta64(int(dt), 's') for dt in time])
        t = ((dates - np.datetime64(date_start)) / np.timedelta64(1, "s"))
        # Loop on inflow sections
        iic = 0
        counter = 0
        for ics in coupling_inflow["ics"].unique():
            if ics == 0:
                continue
            coupling_selection = coupling_inflow[coupling_inflow["ics"] == ics]
            # Load discharge for first hydro_index
            Q = np.zeros(t.size)
            for coupling_index in coupling_selection.index:
                iseg = coupling_selection.loc[coupling_index, "segment"] - 1
                ibase = coupling_selection.loc[coupling_index, "ics"] - 1
                x = coupling_inflow.loc[coupling_index, "x"]
                y = coupling_inflow.loc[coupling_index, "y"]
                Q += data.variables["Q"][:, iic]
                iic += 1
            Q = np.maximum(model.qeps, Q)
            # print(ics, time[0:10], Q[0:10])
            logger.debug("Add inflow %i: (seg=%i, (%.4f,%.4f), Q=[%.3f, %.3f]" % (counter, iseg, x, y, np.min(Q), np.max(Q)))
            model.add_inflow_condition(iseg=iseg, coords=[x, y], t=t, q=Q)
            # plt.plot(t / 86400.0, Q)
            # plt.xlim(0.0, None)
            # plt.title("lat %i, seg:%i, (%.4f,%.4f)" % (counter, iseg, x, y))
            # plt.savefig("plot/direct/debug/qlat%04i.png" % counter)
            # plt.close(plt.gcf())

            counter += 1

        logger.info("Number of inflow points: %i" % iic)
    
    elif "inflow_conditions" in config["model"]:
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
        logger.debug("Number of inflow points: %i" % iic)

    # Set model scheme
    if "numerical_scheme" in config["model"]:
        model.set_scheme(config["model"]["numerical_scheme"])
    else:
        model.set_scheme("preissmann")

    # Set computation timestep
    if "timestep" in config["model"]:
        model.dt = config["model"]["timestep"]
    elif "dt" in config["model"]:
        model.dt = config["model"]["dt"]
    else:
        model.dt = 300
    
    # Set output timestep
    model.dtout = get_from_config(config, section="model", 
                                  key="output-timestep", required=True, 
                                  value_type=float,
                                  deprecated_keys=["output_timestep"])
    # model.dtout = config["model"]["output_timestep"]
    model.output_resampled_cs = True

    output_format = get_from_config(config, section="model", 
                                    key="output-format", required=False, 
                                    value_type=str,
                                    deprecated_keys=["output_format"])
    if output_format is not None:
        if output_format not in ["csv", "netcdf"]:
            raise ValueError("Output format must be either 'csv' or 'netcdf'")
    else:
        output_format = "netcdf"
    model.output_format = output_format
    
    # Set output file
    output_file = get_from_config(config, section="model", 
                                  key="output-file", required=False, 
                                  value_type=str,
                                  deprecated_keys=["output_file"])
    if output_file is not None:
        model.output_file = output_file
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
        
    ## Set gamma for the regularization
    #model.gamma_reg = 0.0
    
    # Set implicit theta for Preissmann scheme
    model.theta_preissmann = get_from_config(config, section="model", 
                                  key="theta-preissmann", required=False, 
                                  value_type=float, default=0.667,
                                  deprecated_keys=["theta_preissmann"])
    # if "theta_preissmann" in config["model"]:
    #     model.theta_preissmann = config["model"]["theta-preissmann"]
    # else:
    #     model.theta_preissmann = 0.667
        
    if "frLPI" in config["model"]:
        model.frLPI = config["model"]["frLPI"]
    else:
        model.frLPI = 2.0
        
    if "checkpoint" in config["model"]:
        model.checkpoint_nt = config["model"]["checkpoint"]
    else:
        model.checkpoint_nt = 5000
            
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
    date_start = get_from_config(config, section="model", 
                                 key="start-date", required=False, 
                                 value_type=validated_date_string, 
                                 deprecated_keys=["date_start"])
    # if "date_start" in config["model"]:
    #     date_start = config["model"]["date_start"]
    # else:
    #     date_start = None

    assim_start_t = 0.0 
    if "assimilation" in config:
        if "date_start" in config["assimilation"]:
            assim_start_date = config["assimilation"]["date_start"]
            assim_start_t = (np.datetime64(assim_start_date) - np.datetime64(date_start)) / np.timedelta64(1, "s")

    # Load metadata file
    if os.path.isfile(os.path.join(static_data_dir, "xs.shp")):
        xs_metadata = gpd.read_file(os.path.join(static_data_dir, "xs.shp"))
        raise NotImplementedError("Finding cross-section index from xs.shp is not implemented yet.")
    elif os.path.isfile(os.path.join(static_data_dir, "xs_nodes.shp")):
        logger.debug("Use xs<->nodes assocation file for observations %s")
        xs_metadata = gpd.read_file(os.path.join(static_data_dir, "xs_nodes.shp"))
    else:
        xs_metadata = None

    # Count number of stations
    nobs = 0
    for iobs, obs_config in enumerate(config["observations"]):
        if "timeseries_file" in obs_config:

            nobs += 1

        elif "provider" in obs_config:

            if obs_config["provider"] == "swot":

                logger.debug("Search for SWOT files")
                logger.debug("- Pattern: %s" % os.path.join(assim_data_dir, "swot", "*_SWOT.nc"))
                candidates_files = glob.glob(os.path.join(assim_data_dir, "swot", "*_SWOT.nc"))
                logger.debug("- Number of SWOT nc files found: %i" % len(candidates_files))

                # Load xs_node file
                # get_from_config(config, "")
                if "nodes-file" in obs_config:
                    xs_nodes_file = obs_config["nodes-file"]
                    xs_nodes_file = auto_filepath(xs_nodes_file, static_data_dir)
                    if not os.path.isfile(xs_nodes_file):
                        logger.error("xs<->nodes file not found: %s" % xs_nodes_file)
                    xs_nodes = gpd.read_file(xs_nodes_file)
                    # print("xs_nodes=")
                    # print(xs_nodes)
                    # choice = input()
                elif xs_metadata is None:
                    xs_nodes_file = "xs_nodes.shp"
                    logger.error("Unable to find xs<->nodes file for provider 'swot'")
                    raise IOError("Unable to find xs<->nodes file for provider 'swot'")
                else:
                    xs_nodes = xs_metadata

                for swot_file in candidates_files:
                    swot_data = nc.Dataset(swot_file, "r")
                    nodes_id = swot_data.groups["node"]["node_id"][:]
                    available_nodes = nodes_id[np.isin(nodes_id, xs_nodes["node_id"].values)]
                    # print("  %i/%i available nodes for file %s" % (available_nodes.size, nodes_id.size, os.path.basename(swot_file)))
                    # print("  ", available_nodes)
                    # choice = input()
                    nobs += available_nodes.size

                logger.debug("- Number of SWOT observations nodes: %i" % nobs)

            elif obs_config["provider"] == "hydroweb":

                logger.debug("Search for hydroweb files")
                candidates_files = glob.glob(os.path.join(assim_data_dir, "hydroweb", "hydroweb_*.txt"))
                logger.debug("- Number of hydroweb files found: %i" % len(candidates_files))
                nobs += len(candidates_files)

            elif obs_config["provider"] == "dahiti-nc":

                dahiti_metadata = gpd.read_file(os.path.join(obs_config["stations-shapefile"]))
                nobs += dahiti_metadata.index.size

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
    total_obs_count = 0
    for iobs, obs_config in enumerate(config["observations"]):
        
        # Load timeseries
        if "timeseries_file" in obs_config:

            tobs, Hobs = load_timeseries(obs_config["timeseries_file"], datestart=date_start)
        
            # Create array of observed H and W, restricted to the assimilation window
            Hobs = Hobs[np.logical_and(tobs >= assim_start_t, tobs <= model.te)]
            tobs = tobs[np.logical_and(tobs >= assim_start_t, tobs <= model.te)]
            HWobs = np.ones((2, Hobs.size)) * -1e+99
            HWobs[0, :] = Hobs
            
            # Setup station
            obs.stations[iobs].setup(model.msh, tobs, HWobs, indices=obs_config["index"])
            ista += 1

        elif "provider" in obs_config:

            if obs_config["provider"] == "swot":

                logger.debug("Load SWOT observations")
                candidates_files = glob.glob(os.path.join(assim_data_dir, "swot", "*_SWOT.nc"))

                # Load xs_node file
                if "nodes-file" in obs_config:
                    xs_nodes_file = obs_config["nodes-file"]
                    xs_nodes_file = auto_filepath(xs_nodes_file, static_data_dir)
                    if not os.path.isfile(xs_nodes_file):
                        logger.error("xs<->nodes file not found: %s" % xs_nodes_file)
                    xs_nodes = gpd.read_file(xs_nodes_file)
                elif xs_metadata is None:
                    xs_nodes_file = "xs_nodes.shp"
                    logger.error("Unable to find xs<->nodes file for provider 'swot'")
                    raise IOError("Unable to find xs<->nodes file for provider 'swot'")
                else:
                    xs_nodes = xs_metadata

                # ista = 0
                # start = 0
                obs_count = 0
                for swot_file in candidates_files:
                    swot_data = nc.Dataset(swot_file, "r")
                    nodes_id = swot_data.groups["node"]["node_id"][:]
                    available_node = np.isin(nodes_id, xs_nodes["node_id"])
                    available_node_id = nodes_id[available_node]
                    sorted_indices = np.argsort(available_node_id)[::-1]
                    available_node_id = available_node_id[sorted_indices]
                    time = swot_data.groups["node"]["time"][sorted_indices, :]
                    # print(time, np.any(np.isnan(time)), time.mask)
                    # choice = input()
                    node_wse = swot_data.groups["node"]["wse"][sorted_indices, :]
                    node_wse[np.isnan(node_wse)] = -1e+99

                    for i in range(0, node_wse.shape[0]):

                        # Filter roughly
                        node_wse
                        node_wse_mean = np.mean(node_wse)
                        # print("[0]", node_wse)
                        node_wse[np.abs((node_wse - node_wse_mean) > 10.0)] = np.nan
                        # print("[1]", node_wse)

                        delta_seconds = (np.datetime64(date_start) - np.datetime64("2000-01-01")) / np.timedelta64(1, "s")
                        # time_node = np.array([np.datetime64("2000-01-01") + np.timedelta64(int(dt), "s") for dt in time[i, :]])
                        # tobs = (time_node - np.datetime64(date_start)) / np.timedelta64(1, "s")
                        tobs = time[i, :] - delta_seconds
                        tobs[time[i, :].mask] = -9999.0
                        valid = np.logical_and(tobs >= 0.0, tobs <= model.te)
                        # valid[41] = False
                        # print("TIME:", np.datetime64(date_start), time_node[41], tobs[41], valid[41])
                        # chcoie = input()
                        # print(time_node - np.datetime64(date_start))
                        # print(valid)
                        # print(time_node)
                        # print(date_start)
                        # if ista < 10:
                            # end = start + np.sum(valid)
                            # print("NODE_ID: %i, ista=%i, wse=[%f, %f], start=%i, end=%i" % (available_node_id[i], ista, np.min(node_wse[i, valid]), np.max(node_wse[i, valid]), start, end))
                            # if ista == 0:
                            #     print(node_wse[i, valid])

                        if np.sum(valid) == 0:
                            tobs = np.ones(1) * 0.0
                            HWobs = np.ones((2, 1)) * -1e+99
                        else:
                            tobs = tobs[valid]
                            HWobs = np.ones((2, np.sum(valid))) * -1e+99
                            # if np.max(node_wse[i, valid]) < 0.0:
                            #     plt.plot(time_node[valid], node_wse[i, valid], "b.")
                            #     plt.show()
                            HWobs[0, :] = node_wse[i, valid]
                        candidates = xs_nodes[xs_nodes["node_id"] == available_node_id[i]]
                        if candidates.index.size == 1:
                            # xs_index = candidates["xs_idx"].values[0]
                            xs_index = candidates.index[0]
                        elif candidates.index.size > 1:
                            # xs_index = candidates["xs_idx"].values[0]
                            xs_index = candidates.index[0]
                        else:
                            raise RuntimeError("%i xs found for node %i" % (candidates.index.size, available_node_id[i]))
                        
                        # Restrict to assimilation window
                        # print(tobs, assim_start_t, model.te)
                        # choice = input()
                        HWobs = HWobs[:, np.logical_and(tobs >= assim_start_t, tobs <= model.te)]
                        tobs = tobs[np.logical_and(tobs >= assim_start_t, tobs <= model.te)]
                        valid_obs = np.sum(HWobs.flatten() > -9999.0)
                        obs_count += valid_obs

                        obs.stations[ista].setup(model.msh, tobs, HWobs, indices=xs_index)
                        obs.stations[ista]._name = str(available_node_id[i])
                        ista += 1

                logger.debug("- Number of SWOT observations: %i" % obs_count)
                total_obs_count += obs_count

            elif obs_config["provider"] == "hydroweb":

                logger.debug("Load hydroweb observations files")
                obs_files = glob.glob(os.path.join(assim_data_dir, "hydroweb", "hydroweb_*.txt"))

                obs_count = 0
                for obs_file in obs_files:

                    logger.debug("- Load hydroweb observations file: %s" % os.path.basename(obs_file))
                    tobs, Hobs, metadata = load_provider_timeseries(obs_file, datestart=date_start)

                    # Restrict to assimilation window
                    Hobs = Hobs[np.logical_and(tobs >= assim_start_t, tobs <= model.te)]
                    tobs = tobs[np.logical_and(tobs >= assim_start_t, tobs <= model.te)]
                    obs_count += tobs.size

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

                logger.debug("- Number of hydroweb observations: %i" % obs_count)
                total_obs_count += obs_count
            
            elif obs_config["provider"] == "dahiti-nc":

                dahiti_metadata = gpd.read_file(os.path.join(obs_config["stations-shapefile"]))

                obs_count = 0
                for index in dahiti_metadata.index:
                    # print(index, dahiti_metadata.loc[index, "dahiti_id"])
                    obs_file = os.path.join(obs_config["data-dir"], "%s_water_level_altimetry.nc" % int(dahiti_metadata.loc[index, "dahiti_id"]))
                    obs_data = nc.Dataset(obs_file, "r")
                    obs_data.set_auto_mask(False)
                    # print(obs_data["datetime"])
                    tobs = np.array([(np.datetime64(date) - np.datetime64(date_start)) / np.timedelta64(1, "s") for date in obs_data["datetime"]])
                    # print(obs_data["wse"][:])
                    egm96 = dahiti_metadata.loc[index, "deltah"]
                    wse = obs_data["wse"][:]
                    wse[~np.isfinite(wse)] = -1e+9

                    Hobs = wse -egm96 - obs_data.geoid
                    # plt.plot(Hobs)

                    # Restrict to assimilation window
                    Hobs = Hobs[np.logical_and(tobs >= assim_start_t, tobs <= model.te)]
                    tobs = tobs[np.logical_and(tobs >= assim_start_t, tobs <= model.te)]
                    # Hobs = Hobs[tobs <= model.te]
                    # tobs = tobs[tobs <= model.te]
                    obs_count += tobs.size

                    HWobs = np.ones((2, Hobs.size)) * -1e+99
                    HWobs[0, :] = Hobs

                    x = dahiti_metadata.loc[index, "longitude"]
                    y = dahiti_metadata.loc[index, "latitude"]
                    logger.debug("  - Station coordinates: %f,%f" % (x, y))
                    obs.stations[ista].setup(model.msh, tobs, HWobs, coords=[x, y])
                    obs.stations[ista]._deltah = -egm96 - obs_data.geoid
                    obs.stations[ista]._egm96 = egm96
                    obs.stations[ista]._geoid = obs_data.geoid
                    obs.stations[ista]._name = str(obs_data.dahiti_id)
                    logger.debug("  - Station index: %i" % obs.stations[ista].ics[0])

                    ista += 1

                logger.debug("- Number of DAHITI observations: %i" % obs_count)
                total_obs_count += obs_count

            elif obs_config["provider"] == "schapi":

                logger.debug("Load schapi observations files")
                if "code_stations" in obs_config:
                    obs_files = [os.path.join(assim_data_dir, "schapi", "schapi_%s.txt" % code) for code in obs_config["code_stations"]]
                else:
                    obs_files = glob.glob(os.path.join(assim_data_dir, "schapi", "schapi_*.txt"))

                obs_count = 0
                for obs_file in obs_files:
                    logger.debug("- Load schapi observations file: %s" % os.path.basename(obs_file))
                    tobs, Hobs, metadata = load_provider_timeseries(obs_file, datestart=date_start)

                    # Restrict to assimilation window
                    Hobs = Hobs[np.logical_and(tobs >= assim_start_t, tobs <= model.te)]
                    tobs = tobs[np.logical_and(tobs >= assim_start_t, tobs <= model.te)]
                    # Hobs = Hobs[tobs <= model.te]
                    # tobs = tobs[tobs <= model.te]
                    obs_count += tobs.size

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

                logger.debug("- Number of SCHAPI observations: %i" % obs_count)
                total_obs_count += obs_count

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

        
    logger.info("Total number of observations: %i" % total_obs_count)
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
                if model.bc[ibc].id.decode("utf-8") == "discharge":
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

        elif item_config["variable"] == "bathy":
            
            # Add bathy
            control.add_bathy_in_control(model, mesh)
            
        else:
            
            raise RuntimeError("Cannot add unknown variable '%s' in control" % item_config["variable"])

        # Setup Background Error Covariance Matrix (BECM) parameters
        if "sigma" in item_config:
            use_becm = True
            becm_sigmas.append(item_config["sigma"])
        else:
            becm_sigmas.append(None)
        if "correlation_length" in item_config:
            use_becm = True
            becm_lengths.append(item_config["correlation_length"])
        else:
            becm_lengths.append(None)

    logger.info("Control size: %i" % control.x.size)

    if use_becm is True:
        logger.info("Use Background Error Correlation Matrix")
        control.x0 = control.x[:]
        # x_prior = control.x0.copy()
        B = control.zero_block_matrix()

        iblock = 0
        for i, item_config in enumerate(config["control"]):

            logger.debug("- Setup BECM block(s) for %s" % item_config["variable"])
            sigma = becm_sigmas[i]
            length = becm_lengths[i]
            if sigma is not None:
                if length is not None:
                    if item_config["variable"][0:3] == "Qin":
                        for ibc in range(0, len(model.bc)):
                            # print(model.bc[ibc].id.decode("utf-8"), model.bc[ibc].id.decode("utf-8") == "discharge")
                            if model.bc[ibc].id.decode("utf-8") == "discharge":
                                sigma_value = sigma
                                if isinstance(sigma, str):
                                    if sigma.strip()[-1] == "%":
                                        sigma_value = float(sigma.strip()[:-1]) * 0.01 * np.mean(model.bc[ibc].ts.y)
                                B.blocks[iblock].m[:, :] = sigma_value**2 * control.temporal_correlation_array(model.bc[ibc].ts.t, length)
                                # plt.imshow(B.blocks[iblock].m[:, :])
                                # plt.show()
                                logger.debug("  - sigma: %f, length=%f" % (sigma_value, length))
                                iblock += 1
                    elif item_config["variable"] in ["K", "alpha", "beta", "KLOB", "KCH", "KROB", "bathy"]:
                        B.blocks[iblock].m[:, :] = sigma**2 * control.spatial_correlation_array(mesh, length)
                        logger.debug("  - sigma: %f, length=%f" % (sigma, length))
                        iblock += 1
                else:
                    if item_config["variable"][0:3] == "Qin":
                        for ibc in range(0, len(model.bc)):
                            if model.bc[ibc].id.decode("utf-8") == "discharge":
                                sigma_value = sigma
                                if isinstance(sigma, str):
                                    if sigma.strip()[-1] == "%":
                                        sigma_value = float(sigma.strip()[:-1]) * 0.01 * np.mean(model.bc[ibc].ts.y)
                                B.blocks[iblock].m[:, :] = sigma_value**2 * np.eye(B.blocks[iblock].m.shape[0])
                                # plt.imshow(B.blocks[iblock].m[:, :])
                                # plt.show()
                                logger.debug("  - sigma: %f" % sigma_value)
                                iblock += 1
                    elif item_config["variable"] in ["K", "alpha", "beta", "KLOB", "KCH", "KROB", "bathy"]:
                        B.blocks[iblock].m[:, :] = sigma**2 * np.eye(B.blocks[iblock].m.shape[0])
                        iblock += 1
            elif length is not None:
                raise ValueError("Cannot use correlation length without sigma")

        control.set_prior_error_covariance_matrixblock(B)
        logger.debug("- Compute Choleski decomposition")
        control.x[:] = 0.0
    logger.info("")
    
    return control

def write_netcdf_results(model, config):

    logger = logging.getLogger("dassflow1d")
    logging.info("Export results to NetCDF file")

    # Compute dates array
    date_start = get_from_config(config, section="model", 
                                 key="start-date", required=True, 
                                 value_type=validated_date_string, 
                                 deprecated_keys=["date_start"])
    date_start = np.datetime64(date_start)
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
                indices.append(ics)

    # Remove outliers in discharge timeseries
    dt = (dates - np.datetime64("2000-01-01 00:00:00")) / np.timedelta64(1, "s")
    Qout = model.res.q[indices, :].T
    for ix in range(0, nx):
        for it in range(0, nt):
            Qwindow = Qout[max(0, it-5):it+5, ix]
            Qmean = np.mean(Qwindow)
            if Qout[it, ix] > 5.0 * Qmean:
                Qout[it, ix] = Qmean
            # valid = Qout[:, ix]/Qmean <= 5.0
            # unvalid = Qout[:, ix]/Qmean > 5.0
            # if np.sum(valid) > 0:
            #     Qout[unvalid, ix] = np.interp(dt[unvalid], dt[valid], Qout[valid, ix])

    # Create netCDF output file
    logger.info("- Create NetCDF output file: %s" % model.output_nc_file)
    dataset = nc.Dataset(model.output_nc_file, "w", format="NETCDF4")
    dataset.createDimension("nx", nx)
    dataset.createDimension("nt", nt)
    time = dataset.createVariable("time", "f4", ("nt",))
    time.setncattr("units", "seconds since 2000-01-01 00:00:00")
    time[:] = dt
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
    Q[:, :] = Qout
    dataset.close()





def run_direct(config, display_internal_counters=False, obs=True):

    logger = logging.getLogger("dassflow1d")
  
    # Create model
    model = create_model(config)

    if obs is True:
        obs = load_observations(config, model)
    else:
        obs = None
    
    # Run direct
    logger.info("=" * 80)
    logger.info(" SIMULATION")
    logger.info("=" * 80)
    
    if obs is None:
        # Run direct model
        model.run_unsteady()
    else:
        # Run calc cost with a dummy control to generate estimation
        control = m_control.Control()
        cost = dassflow1d.calc_cost(model, control, obs)
        logger.info("status=%i" % model.status)
        logger.info("Cost=%i" % cost)


    if model.output_format == "netcdf":
        write_netcdf_results(model, config)

    if display_internal_counters:
        logger.debug("- Number of overbank flow occurences: %i" % model.internal_counters[0])
        logger.debug("  - Left overbank flow occurences : %i" % model.internal_counters[1])
        logger.debug("  - Right overbank flow occurences: %i" % model.internal_counters[2])
        logger.debug("- Number of overflow (flow above highest elevation) occurences: %i" % model.internal_counters[3])
    logger.info("")

    # for iseg in range(0, model.msh.nseg):

    #     x = model.msh.get_segment_field(iseg, "x")
    #     bathy = model.msh.get_segment_field(iseg, "bathy")
    #     z = model.get_segment_results(iseg, 0, "z")
    #     plt.plot(x, bathy, "k-")
    #     plt.plot(x, z, "b-")
    #     plt.show()
    if "plots" in config:
        plots_config = config["plots"]
        if "results-prefix" in plots_config:
            plot_results_prefix = plots_config["results-prefix"]
            plot_results = True
            if "plot-results" in plots_config:
                plot_results = plots_config["plot-results"]
            if plot_results:
                for iout in range(0, model.res.t.size, 8):
                    date = (datetime(2019, 1, 1) + timedelta(seconds=model.res.t[iout])).strftime("%Y%m%dT%H%M")
                    print("Plot results for %s (%03i %%)" % (date, (iout+1) * 100 // model.res.t.size))
                    for iseg in range(0, model.msh.nseg):
                        # Retrieve field 'x' (distance from outlet) for segment 0
                        x = model.msh.get_segment_field(iseg=iseg, field="x")
                        # Retrieve field 'bathy' (bathymetry) for segment 0
                        bathy = model.msh.get_segment_field(iseg=iseg, field="bathy")
                        # Retrieve results field 'z' (water surface elevation) for segment 0 at output index 0
                        z = model.get_segment_results(iseg=iseg, iout=iout, field="z")
                        # Retrieve results field 'q' (discharge) for segment 0 at output index 0
                        q = model.get_segment_results(iseg=iseg, iout=iout, field="q")
                        # Retrieve results field 'w' (width) for segment 0 at output index 0
                        w = model.get_segment_results(iseg=iseg, iout=iout, field="w")
                        # Create and show plot
                        plot_file = "%s_seg%03i_out%08i.png" % (plot_results_prefix, iseg+1, iout)
                        plot_BZQW(x, z, q, w, time_label="%s (reach %i)" % (date, iseg), bathy=bathy, outfile=plot_file)


    if obs is not None:
        
        for iseg in range(0, model.msh.nseg):

            err_seg = np.zeros(0)

            for ista in range(0, len(obs.stations)):    

                ics = obs.stations[ista].ics[0]
                if ics >= model.msh.seg[iseg].first_cs and ics <= model.msh.seg[iseg].last_cs:
                    start = obs.stations[ista].offset
                    end = start + len(obs.stations[ista].t)
                    obs_sta = obs.obs[0, start:end]
                    est_sta = obs.est[0, start:end]
                    valid = np.logical_and(est_sta > -9999.0, obs_sta > -9999.0)
                    obs_sta = obs_sta[valid]
                    est_sta = est_sta[valid]
                    err_sta = est_sta - obs_sta
                    err_seg = np.concatenate((err_seg, err_sta))

            if err_seg.size > 0:

                y, x, _ = plt.hist(err_seg, 100, density=True, cumulative=True, histtype="step", alpha=0.0)
                y = np.insert(y, 0, 0.0)
                
                plt.plot(x, y, "k-", label="prior")
                plt.axvline(x[0], color="k", ls="--")
                plt.axvline(x[-1], color="k", ls="--")
                plt.xlabel(r"$H - H_obs$")
                plt.ylabel("CDF")
                plt.title("segment %i" % iseg)
                plt.tight_layout()
                plt.show()

        err = np.zeros(0)
        x_sta = []
        Hm_sta = []
        for ista in range(0, len(obs.stations)):
            start = obs.stations[ista].offset
            end = start + len(obs.stations[ista].t)
            obs_sta = obs.obs[0, start:end]
            est_sta = obs.est[0, start:end]
            valid = np.logical_and(est_sta > -9999.0, obs_sta > -9999.0)
            obs_sta = obs_sta[valid]
            est_sta = est_sta[valid]
            Hm_sta.append(np.mean(obs_sta))
            err_sta = est_sta - obs_sta
            err = np.concatenate((err, err_sta))
            #err = np.concatenate((err, obs.est[0, start:end] - obs.obs[0, start:end]))
            # if ista < 100:
            #     print(ista, np.min(err_sta), np.max(err_sta))
            #     if ista == 0:
            #         print(obs.obs[0, start:end])
            #         print(obs.est[0, start:end])

            # if obs.stations[ista]._name == "12073":
            bathy = model.msh.cs[obs.stations[ista].ics[0]].bathy
            plt.plot(obs_sta, "r.")
            plt.plot(est_sta, "b+")
            plt.axhline(bathy, c="k", ls="--")
            plt.title("Station %i (ics:%s), name=%s" % (ista, str(obs.stations[ista].ics), obs.stations[ista]._name))
            print(model.msh.cs[obs.stations[ista].ics[0]].ibase)
            plt.show()
        Hm_sta = np.array(Hm_sta)
        Hm_sta[Hm_sta <= -9999] = np.nan
        plt.plot(Hm_sta, "r.")
        plt.show()
        Hm_sta_valid = Hm_sta[Hm_sta >= 0.0]
        print("Hm_sta_valid: count=%i" % Hm_sta_valid.size)
        plt.plot(Hm_sta_valid, "r.")
        plt.show()
        plt.hist(err, bins=100)
        plt.show()

        #     if err.size > 0:
        #         figfile = os.path.join("plot", "err_prior_segment%i.png" % iseg)
        #         fig = plt.figure()
        #         y, x, _ = plt.hist(err, 100, density=True, cumulative=True, histtype="step", alpha=0.0)
        #         y = np.insert(y, 0, 0.0)
        #         plt.plot(x, y, "k-", label="prior")
        #         plt.axvline(x[0], color="k", ls="--")
        #         plt.axvline(x[-1], color="k", ls="--")
        #         plt.xlabel(r"$H - H_obs$")
        #         plt.ylabel("CDF")
        #         plt.legend()
        #         plt.tight_layout()
        #         plt.savefig(figfile)
        #         plt.close(fig)

        # return err_prior


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
    
class MinimizeCallback:

    def __init__(self, config, model, control=None, obs=None, output_file=None):
        self.config = config
        self.model = model
        self.control = control
        self.obs = obs
        x_prior = np.zeros(control.x.size)
        control.get_model_control(model, x_prior)
        self.data = {
            "cost": np.zeros(0),
            "grad": np.zeros((0, control.x.size)),
            "ctrl-x": np.zeros((0, control.x.size)),
            "ctrl-x0": x_prior}
        self.iteration = 0
        self._plot_obs = False
        self._output_file = output_file

    def __call__(self, x):
        self.iteration += 1
        x_iter = np.zeros(x.size)
        self.control.get_model_control(self.model, x_iter)
        cost = self.model.__last_cost__
        grad = self.model.__last_grad__
        logger.info("Iteration %03i, cost= %12.5e, [proj g]= %12.5e" % (self.iteration, cost, np.linalg.norm(grad)))
        self.data["cost"] = np.concatenate((self.data["cost"], [cost]))
        self.data["grad"] = np.concatenate((self.data["grad"], grad.reshape((1,-1))))
        self.data["ctrl-x"] = np.concatenate((self.data["ctrl-x"], x_iter.reshape((1,-1))))

        # Save optimization results
        # out_file = None
        # if "assimilation" in self.config:
        #     assimilation_config = self.config["assimilation"]
        #     if "output-file" in assimilation_config:
        #         out_file = assimilation_config["output-file"]
        # elif "optim" in self.config:
        #     optim_config = self.config["optim"]
        #     if "output-file" in optim_config:
        #         out_file = optim_config["output-file"]
        # if out_file is not None:
        if self._output_file is not None:
            out = nc.Dataset(self._output_file, "w")
            out.createDimension("nchars", 6)
            out.createDimension("nitems", len(self.control.items))
            out.createDimension("nx", self.control.x.size)
            out.createDimension("niter", self.iteration)
            items = out.createVariable("items", 'S1', ("nitems", "nchars"))
            for iitem, item in enumerate(self.control.items):
                for ichar in range(0, 6):
                    if ichar < len(item.id):
                        items[iitem, ichar] = item.id[ichar:ichar+1]
                    else:
                        items[iitem, ichar] = '\0'
            items_offset = out.createVariable("items_offset", 'i4', ("nitems",))
            for iitem, item in enumerate(self.control.items):
                items_offset[iitem] = item.offset
            items_size = out.createVariable("items_size", 'i4', ("nitems",))
            for iitem, item in enumerate(self.control.items):
                items_size[iitem] = item.nx
            out_x0 = out.createVariable("x0", "f4", ("nx"))
            out_x0[:] = self.data["ctrl-x0"][:]
            out_x = out.createVariable("x", "f4", ("niter", "nx"))
            out_x[:, :] = self.data["ctrl-x"][:, :]
            out_cost = out.createVariable("cost", "f4", ("niter"))
            out_cost[:] = self.data["cost"][:]
            out_grad = out.createVariable("grad", "f4", ("niter", "nx"))
            out_grad[:, :] = self.data["grad"][:, :]

            out.close()

        if self._plot_obs is True:
            if os.path.isdir("plot/optim"):
                self.plot_obs_cdf()
                self.plot_obs_ts()

    def plot_obs_cdf(self):

        mesh = self.model.msh
        for iseg in range(0, mesh.nseg):
            print("SEG:", iseg)

            # Select stations
            err_h_segment = np.zeros(0)
            for ista in range(0, len(self.obs.stations)):
                # print("sta %05i, ista=%s seg=[%i - %i]" % (ista, str(self.obs.stations[ista].ics), mesh.seg[iseg].first_cs, mesh.seg[iseg].last_cs))
                if np.any(np.logical_and(self.obs.stations[ista].ics >= mesh.seg[iseg].first_cs,
                                         self.obs.stations[ista].ics <= mesh.seg[iseg].last_cs)):
                    start = self.obs.stations[ista].offset
                    end = start + len(self.obs.stations[ista].t)
                    obs_sta = self.obs.obs[0, start:end]
                    est_sta = self.obs.est[0, start:end]
                    valid = np.logical_and(est_sta > -9999.0, obs_sta > -9999.0)
                    # obs_sta[obs_sta < -9999.0] = np.nan
                    # est_sta[est_sta < -9999.0] = np.nan
                    # valid = np.logical_and(station_est > -1e+9, station_obs > -1e+9)
                    # print("valid: %i (obs:%i, est:%i)" % (np.sum(valid), np.sum(obs_sta > -9999.0), np.sum(est_sta > -9999.0)))
                    # choice = input()
                    if np.sum(valid) > 0:
                        err_h_station = est_sta[valid] - obs_sta[valid]
                        err_h_segment = np.concatenate((err_h_segment, err_h_station))

            # print(" => err_h_segment.size=%i" % err_h_segment.size)
            if err_h_segment.size > 0:
                # fig = plt.figure()
                y, x, _ = plt.hist(err_h_segment, 100, density=True, cumulative=True, histtype="step", alpha=0.0)
                y = np.insert(y, 0, 0.0)
                plt.plot(x, y, "b-", label="prior")
                plt.axvline(x[0], color="b", ls="--")
                plt.axvline(x[-1], color="b", ls="--")
                plt.xlabel("H - Hobs")
                plt.ylabel("CDF")
                plt.tight_layout()
                figfile = "plot/optim/iter%04i_segment%04i_errH_CDF.png" % (self.iteration, iseg)
                # print(figfile)
                try:
                    plt.savefig(figfile)
                except:
                    plt.show()
                    raise
                plt.close(plt.gcf())
                # plt.show()

    def plot_obs_ts(self):

        mesh = self.model.msh


        fig = None
        panel_count = 0
        for ista in range(0, len(self.obs.stations)):

            if ista%16 == 0:
                if fig is not None:
                    figfile = "plot/optim/iter%04i_stations_panel%04i.png" % (self.iteration, panel_count)
                    plt.savefig(figfile)
                    plt.close(fig)
                    fig = None
                    panel_count += 1
                fig = plt.figure()
                plotgrid = gs.GridSpec(4, 4)

            igs = ista%16
            ax = fig.add_subplot(plotgrid[igs])
            start = self.obs.stations[ista].offset
            end = start + len(self.obs.stations[ista].t)
            obs_sta = self.obs.obs[0, start:end]
            est_sta = self.obs.est[0, start:end]
            valid = np.logical_and(est_sta > -9999.0, obs_sta > -9999.0)
            obs_sta[obs_sta < -9999.0] = np.nan
            est_sta[est_sta < -9999.0] = np.nan
            bathy = np.mean([mesh.cs[ics].bathy for ics in self.obs.stations[ista].ics])
            ax.plot(obs_sta, "r.", label="obs")
            ax.plot(est_sta, "b+", label="est")
            if igs%4 == 0:
                ax.set_ylabel("H (m)")
            ax.axhline(bathy, c="k", ls="--")


        if fig is not None:
            # ax.suptitle("Station %i (ics:%s), name=%s" % (ista, str(obs.stations[ista].ics), obs.stations[ista]._name))
            figfile = "plot/optim/iter%04i_stations_panel%04i.png" % (self.iteration, panel_count)
            plt.savefig(figfile)
            plt.close(fig)
            fig = None


def run_optim(config, max_iterations=200, feps=1e-3):

    logger = logging.getLogger("dassflow1d")

    # Create model
    model = create_model(config)
  
    # Load observations
    obs = load_observations(config, model)
    if len(obs.stations) == 0:
        raise RuntimeError("No observations stations")
  
    # Setup control
    if "control" not in config:
        raise RuntimeError("Section 'control' not found in config")
    if len(config["control"]) == 0:
        raise RuntimeError("No control items in config")
    control = load_control(config, model)

    # Delete existing output files
    # TODO implement restart
    optim_output_file = get_from_config(config, section="assimilation", 
                                  key="output-file", required=False, 
                                  value_type=str,
                                  deprecated_keys=["output_file"])
    if optim_output_file is not None:
        if os.path.isfile(optim_output_file):
            os.remove(optim_output_file)
    output_file = get_from_config(config, section="model", 
                                  key="output-file", required=False, 
                                  value_type=str,
                                  deprecated_keys=["output_file"])
    if output_file is not None:
        if os.path.isfile(output_file):
            os.remove(output_file)
    
    # Run optimisation
    logger.info("=" * 80)
    logger.info(" OPTIMISATION")
    logger.info("=" * 80)
    
    # Run adjoint model to compute cost and gradient
    dtout0 = model.dtout
    model.dtout = -1
    model.disable_stdout = True

    iter_costs = np.zeros(0)
    iter_grads = np.zeros((0, control.x.size))
    iter_control = np.zeros((0, control.x.size))

    # def callback_minimize(x):
    #     global iter_costs
    #     x_iter = np.zeros(x.size)
    #     control.get_model_control(model, x_iter)
    #     cost = model.__last_cost__
    #     grad = model.__last_grad__
    #     iter_costs = np.concatenate((iter_costs, [cost]))
    #     iter_grads = np.concatenate((iter_grads, grad.shape((1,-1))))
    #     iter_control = np.concatenate((iter_grads, x_iter.shape((1,-1))))

    #     # Save optimization results
    #     if "optim" in config:
    #         if "output-dir" in config["optim"]:
    #             out_file = os.path.join(config["optim"]["output-dir"], "optim.nc")
    #             out = nc.Dataset(out_file, "w")
    #             out.createDimension("nchars", 6)
    #             out.createDimension("nitems", len(control.items))
    #             out.createDimension("nx", control.x.size)
    #             out.createDimension("niter", iter_costs.size)
    #             items = out.createVariable("items", 'S1', ("nitems", "nchars"))
    #             items[:] = nc.stringtochar(np.array([item.id for item in control.items]))

    #             out_x = out.createVariable("x", "f4", ("niter", "nx"))
    #             out_x[:, :] = iter_control[:, :]
    #             out_cost = out.createVariable("cost", "f4", ("niter"))
    #             out_cost[:] = iter_costs[:]
    #             out_grad = out.createVariable("grad", "f4", ("niter", "nx"))
    #             out_grad[:, :] = iter_grads[:, :]


    #     logger.info("Iteration, cost= %12.5e, [proj g]= %12.5e" % (model.__last_cost__, np.linalg.norm(model.__last_grad__)))

    callback_minimize = MinimizeCallback(config, model, control=control, obs=obs, output_file=optim_output_file)

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
    parser.add_argument("-root-dir", dest="rootdir", default=None, 
                        help="Path to the root directory")
    parser.add_argument("--display-internal-counters", dest="display_internal_counters", action="store_true", 
                        help="Display internal counters values")
    parser.add_argument("--debug", action="store_true", 
                        help="Enable debug log")
    args = parser.parse_args()

    # Create Logger
    logger = create_logger(args.debug)
    
    # Load configuration
    config = load_configuration(args.config_file)

    if args.rootdir is not None:
        __rootdir__ = args.rootdir
    
    # Process run
    if args.run_type == "direct":

        # Direct run
        withobs = "observations" in config
        run_direct(config, display_internal_counters=args.display_internal_counters, obs=withobs)

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
