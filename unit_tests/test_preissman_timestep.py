import numpy as np
import os
import matplotlib.pyplot as plt

import dassflow1d
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
from dassflow1d.post.results import plot_BZQ


#def test_preissmann_timestep_01(plot=False):
    #""" Test preissmann timestep on mesh mesh_case_forward_channel.geo
    #"""

    ## Read mesh
    #mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    #mesh.set_uniform_strickler_parameters([20.0, 0.0])
    
    ## Create model and unknowns
    #model = m_sw_mono.Model(mesh)
    
    #model.bc[0].id = "discharge"
    #model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    ##model.bc[1].id = "normal_depth"
    #model.bc[1].id = "elevation"
    #model.bc[1].set_timeseries(t=[0.0, 7200.0], y=[1.0, 1.0])
    #model.set_scheme("preissmann")

    ## Standard step
    #model.standard_step(1e-4, 1000)
    
    ## Preissmann timestep
    #model.dt = 60
    #model.theta_preissmann = 0.667
    
    #status = 0
    #while model.tc < 7200.0:
      #model.tc += model.dt
      #model.preissmann_timestep(model.msh, model.imp, model.dof, status)
      ##print("Q=", dof.q)
      ##print("h=", dof.h)
      ##choice = input()
    
    #h_valid = [1.08276526, 1.10308899, 1.12279015, 1.14164275, 1.15930825,
               #1.17527517, 1.18875893, 1.19852757, 1.20257567, 1.19745117,
               #1.17664581, 1.12574787, 1.        ]
    #q_valid = [1.        , 1.02416142, 1.048116  , 1.07182473, 1.09523228,
               #1.11826061, 1.14079926, 1.16268915, 1.18369337, 1.20343709,
               #1.22126054, 1.2357418 , 1.24189833]
    
    #if plot:
      #plt.plot(model.dof.h[2:-2])
      #plt.plot(mesh.bathy[2:-2])
      #plt.show()
    ##print(dof.h[2:-2])
    ##print(h_valid)
    #assert(np.allclose(model.dof.h[2:-2], h_valid, atol=1e-3))
    ##print(dof.q[2:-2])
    ##print(q_valid)
    #assert(np.allclose(model.dof.q[2:-2], q_valid, atol=1e-3))


def test_preissmann_timestep_02(plot=False):
    """ Test preissmann timestep on mesh mesh_case_forward_channel.geo
    """

    # Read mesh
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    
    mesh.set_uniform_strickler_parameters([20.0, 0.0])
    
    # Create model
    model = m_sw_mono.Model(mesh)
    model.bc[0].id = "discharge"
    model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    model.bc[1].id = "normal_depth"
    model.set_scheme("preissmann")

    # Standard step
    model.tc = 0.0
    model.standard_step(1e-4, 1000)
    #print("Q=", model.dof.q)
    #print("h=", model.dof.h)
    assert(np.allclose(model.dof.q, 1.0))
    assert(np.allclose(model.dof.h[2:-2], 1.03270875))
    
    # Preissmann timestep
    model.dt = 60
    model.theta_preissmann = 0.667
    
    status = 0
    while model.tc < 7200.0:
      model.tc += model.dt
      model.preissmann_timestep(model.msh, model.imp, model.dof, status)
      #print("Q=", dof.q)
      #print("h=", dof.h)
      #choice = input()
    
    h_valid = [1.0841429 , 1.10533756, 1.1263077 , 1.14700479, 1.16735269,
               1.18723369, 1.2064665 , 1.22477122, 1.24171196, 1.25659997,
               1.26832156, 1.27501142, 1.27336724]
    q_valid = [1.        , 1.02437825, 1.04864936, 1.0728104 , 1.09685548,
               1.12077529, 1.14455677, 1.16818314, 1.19163475, 1.21489208,
               1.23794308, 1.26080101, 1.28354836]
    print("Q=", model.dof.q[2:-2])
    print("h=", model.dof.h[2:-2])
    if plot:
      x0 = mesh.get_segment_field(0, "x")
      b0 = mesh.get_segment_field(0, "bathy")
      z0 = model.dof.get_segment_field(mesh, 0, "z")
      q0 = model.dof.get_segment_field(mesh, 0, "q")
      plot_BZQ(x0, z0, q0, time_label="final", bathy=b0)
      
      plt.show()
    assert(np.allclose(model.dof.h[2:-2], h_valid, atol=1e-3))
    assert(np.allclose(model.dof.q[2:-2], q_valid, atol=1e-3))


# To run tests without pytest (debug)
if __name__ == "__main__":
    #test_preissmann_timestep_01()
    #test_preissmann_timestep_02(plot=True)
    pass
