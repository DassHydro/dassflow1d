import numpy as np
import os

import dassflow1d
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_obs as m_obs


# TODO: activate again after #TASK 5557
#def test_strickler_fields_k_constant():

    ## Read mesh
    #mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    #mesh.set_strickler_type("constant")
    #mesh.set_strickler_fields_linear([0.0, 300.0], [[20.0, 25.0]])
    
    #mesh.apply_strickler_fields()
    #ics = 0
    #for cs in mesh.cs:
      ##print("CS(%i):x=%f, K=%f, %f" % (ics, cs.x, cs.strickler_params[0], cs.strickler_params[1]))
      #if ics < 2:
        #assert(cs.strickler_params[0] == 25.0)
        #assert(cs.strickler_params[1] == 0.0)
      #elif ics > 14:
        #assert(cs.strickler_params[0] == 20.0)
        #assert(cs.strickler_params[1] == 0.0)
      #else:
        #assert(np.abs(cs.strickler_params[0] - (25.0 - (ics-2) * 5.0 / 12)) < 1e-6)
        #assert(cs.strickler_params[1] == 0.0)
      #ics += 1


# TODO: activate again after #TASK 5557
#def test_strickler_fields_k_powerlaw_h():

    ## Read mesh
    #mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    #mesh.set_strickler_type("powerlaw_h")
    #mesh.set_strickler_fields_linear([0.0, 300.0], [[20.0, 25.0], [0.2, -0.2]])
    
    #mesh.apply_strickler_fields()
    #ics = 0
    #for cs in mesh.cs:
      #if ics < 2:
        #assert(cs.strickler_params[0] == 25.0)
        #assert(cs.strickler_params[1] == -0.2)
      #elif ics > 14:
        #assert(cs.strickler_params[0] == 20.0)
        #assert(cs.strickler_params[1] == 0.2)
      #else:
        #assert(np.abs(cs.strickler_params[0] - (25.0 - (ics-2) * 5.0 / 12)) < 1e-6)
        #assert(np.abs(cs.strickler_params[1] - (-0.2 + (ics-2) * 0.4 / 12)) < 1e-6)
      #ics += 1


#def test_strickler_fields_01():


    ## Read mesh
    #mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    #mesh.set_uniform_strickler_parameters([20.0, 0.0])
    
    ## Create model and unknowns
    #model = m_sw_mono.Model(mesh)
    #dof = m_sw_mono.Unknowns(mesh)
    #model.bc[0].id = "discharge"
    #model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    #model.bc[1].id = "normal_depth"
    #imp = m_sw_mono.Implicitmatrix(model, mesh)
    #model.ts = 0
    #model.te = 7200
    #model.dt = 60
    #model.theta_preissmann = 0.667

    ## Set inflow
    #dof.q[0:2] = 1.0
    #print(dof.q)

    ## Create observations
    #obs = m_obs.all_observed([0.0,60.0,120.0, 7200.0], mesh)
    
    ## Time loop
    #model.time_loop(mesh, imp, dof, obs)
    #print("GENOBS:")
    #print(dof.h)
    #print(dof.q)
    ##obsdata = m_obs.Obsstation([1,2,3], [0.0,60.0,120.0], [[1.0, 1.0,1.0],[0.1,0.1,0.1]])
    #print(obs.est)
    
    #obs.obs[:, :] = obs.est[:, :]
    #model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[2.0, 2.0, 4.0, 2.0, 1.5])
    #dof.q[0:2] = 1.0
    #model.time_loop(mesh, imp, dof, obs)
    #print("TEST:")
    #print(dof.h)
    #print(dof.q)
    #print(obs.obs)
    #print(obs.est)
    ##print(obs.stations[:].est)
    ##hobs = obs.get_all_observed_heights()
    ##hest = obs.get_all_estimated_heights()
    
    #print ("cost=", np.linalg.norm((obs.obs[0] - obs.est[0])**2))
    
    
    ##print(dof.q)
    
    
    ##print(obs)
    ##print(obs.stations[0])
    ##unk = m_sw_mono.Unknowns(mesh)
    ##assert(unk.a.size == 2)
    ##assert(unk.q.size == 2)
    ##assert(unk.h.size == 2)
    ##assert(unk.sg.size == 2)
    ##assert(unk.sf.size == 2)


# To run tests without pytest (debug)
if __name__ == "__main__":
    test_strickler_fields_k_constant()
    test_strickler_fields_k_powerlaw_h()
