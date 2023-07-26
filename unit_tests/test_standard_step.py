import numpy as np
import os
import matplotlib.pyplot as plt

import dassflow1d
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono


def test_standard_step_01():
    """ Test standard step on a mesh with 2 cross-sections and a single segment
    """

    # Create mesh
    mesh = m_mesh.Mesh(14, 1)
    for i in range(0, 14):
      z = (13 - i) * 0.01
      mesh.setup_crosssection(i+1, 1)
      mesh.cs[i].x = i * 10.0
      mesh.cs[i].level_heights[:] = z + 1.0
      mesh.cs[i].level_widths[:] = 100.0
      mesh.cs[i].bathy = z
      mesh.cs[i].deltademi = 10.0
    mesh.seg[0].first_cs = 3
    mesh.seg[0].last_cs = 12
    mesh.seg[0].ds_seg = 0
    mesh.seg[0].us_bc = 1
    mesh.seg[0].ds_bc = 2
    mesh.set_uniform_strickler_parameters([20.0, 0.0])
    mesh.update_geometries()
    
    # Create model
    model = m_sw_mono.Model(mesh)
    
    # Set boundary conditions
    model.bc[0].id = "discharge"
    model.bc[0].set_timeseries(t=[0.0, 100.0], y=[1.0, 1.0])
    #model.bc[1].id = "elevation"
    #model.bc[1].ts = m_sw_mono.Timeseries(0, 0.2)
    model.bc[1].id = "normal_depth"
    #model.bc[1].ts = m_sw_mono.Timeseries(0, 0.2)

    # Standard step
    model.standard_step(1e-4, 1000)
    h_valid = [2.59200173, 2.60200163, 2.61200152, 2.62200142, 2.63200132,
               2.64200122, 2.65200112, 2.66200102, 2.67200093, 2.68200083]
    assert(np.allclose(model.dof.q, 1.0))
    #print(model.dof.h[2:-2])
    assert(np.allclose(model.dof.h[2:-2], h_valid))


def test_standard_step_02():
    """ Test standard step on mesh mesh_case_forward_channel.geo
    """

    # Read mesh
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    
    mesh.set_uniform_strickler_parameters([20.0, 0.0])
    
    # Create model and unknowns
    model = m_sw_mono.Model(mesh)
    
    # Set boundary conditions
    model.bc[0].id = "discharge"
    model.bc[0].set_timeseries(t=[0.0, 100.0], y=[1.0, 1.0])
    model.bc[1].id = "normal_depth"
    
    # Standard step
    model.standard_step(1e-4, 1000)
    assert(np.allclose(model.dof.q, 1.0))
    assert(np.allclose(model.dof.h[2:-2], 1.03270875))


# To run tests without pytest (debug)
if __name__ == "__main__":
    test_standard_step_01()
