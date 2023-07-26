import numpy as np
import os

import dassflow1d
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_obs as m_obs


def test_time_loop_01():

    # Read mesh
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    mesh.set_uniform_strickler_parameters([20.0, 0.0])
    
    # Create model and unknowns
    model = m_sw_mono.Model(mesh)
    model.bc[0].id = "discharge"
    model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    model.bc[1].id = "normal_depth"
    model.ts = 0
    model.te = 7200
    model.dt = 60
    model.dtout = -1
    model.theta_preissmann = 0.667
    model.set_scheme("preissmann")
    
    # Time loop
    model.time_loop()
    
    h_target = [1.08414187, 1.1053367,  1.1263069,  1.14700404, 1.16735198,
                1.18723301, 1.20646584, 1.22477057, 1.24171132, 1.25659935,
                1.268321,   1.27501107, 1.27336751]
    q_target = [1.,         1.0243783,  1.04864944, 1.07281051, 1.09685562, 
                1.12077545, 1.14455696, 1.16818335, 1.19163499, 1.21489234,
                1.23794338, 1.26080133, 1.28354868]


# To run tests without pytest (debug)
if __name__ == "__main__":
    test_time_loop_01()
