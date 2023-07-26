import numpy as np
import os

import dassflow1d
import dassflow1d.m_control as m_control
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_obs as m_obs


def test_all_observed():

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
    #model.set_scheme("implicit_diffusive_wave")
    # Create observations
    obs = m_obs.all_observed([0.0,60.0,120.0, 7200.0], mesh)
    
    # Time loop
    model.generate_observations(obs)
    
    # Setup control
    control = m_control.Control()
    control.add_bc_in_control(model, 0)
    control.x[0:5] = [2.0, 2.0, 4.0, 2.0, 1.5]
    
    cost = dassflow1d.calc_cost(model, control, obs)
    print('cost',cost, np.allclose(cost, 27.451394309121454))
    #choice = input()
    assert(np.allclose(cost, 27.451394309121454))


# To run tests without pytest (debug)
if __name__ == "__main__":
    test_all_observed()
