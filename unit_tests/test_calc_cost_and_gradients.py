import numpy as np
import os
import sys

import dassflow1d
import dassflow1d.m_control as m_control
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_obs as m_obs


def test_all_observed():


    # Read mesh
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    mesh.set_strickler_type("powerlaw_h")
    mesh.set_uniform_strickler_parameters([20.0, 0.0])
    
    # Create model and unknowns
    model = m_sw_mono.Model(mesh)
    model.bc[0].id = "discharge"
    model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    model.bc[1].id = "normal_depth"
    #imp = m_sw_mono.Implicitmatrix(model, mesh)
    model.ts = 0
    model.te = 7200
    model.dt = 60
    model.dtout = -1
    model.theta_preissmann = 0.667
    model.set_scheme("preissmann")
    #model.set_scheme("implicit_diffusive_wave")

    # Generate observations
    obs = m_obs.all_observed([0.0,60.0,120.0, 7200.0], mesh)
    model.generate_observations(obs)
    
    # Setup control
    control = m_control.Control()
    control.add_bc_in_control(model, 0)
    control.x[0:5] = [2.0, 2.0, 4.0, 2.0, 1.5]
    
    #model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[2.0, 2.0, 4.0, 2.0, 1.5])
    cost, grad = dassflow1d.calc_cost_and_gradients(model, control, obs)
    cost_target = 27.451394309121454
    #grad_target = np.array([5.17981200e+01, 4.14673071e-01, 7.00174593e-04, 5.71318712e-01, 5.73441795e+00])
    grad_target = np.array([5.17981248e+01, 4.14673721e-01, 7.00159669e-04, 5.71388443e-01, 5.73614003e+00])
    print("1:", cost)
    print("2:", grad)
    assert(np.allclose(cost, cost_target))
    assert(np.allclose(grad, grad_target))


def test_all_observed_change_of_variable():


    # Read mesh
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    mesh.set_strickler_type("powerlaw_h")
    mesh.set_uniform_strickler_parameters([20.0, 0.0])
    
    # Create model and unknowns
    model = m_sw_mono.Model(mesh)
    model.bc[0].id = "discharge"
    model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    model.bc[1].id = "normal_depth"
    #imp = m_sw_mono.Implicitmatrix(model, mesh)
    model.ts = 0
    model.te = 7200
    model.dt = 60
    model.dtout = -1
    model.theta_preissmann = 0.667
    #model.set_scheme("implicit_diffusive_wave")
    model.set_scheme("preissmann")
 
    # Generate observations
    obs = m_obs.all_observed([0.0,60.0,120.0, 7200.0], mesh)
    model.generate_observations(obs)
    
    # Setup control
    control = m_control.Control()
    control.add_bc_in_control(model, 0)
    #control.x0[0:5] = [2.0, 2.0, 4.0, 2.0, 1.5]
    control.set_prior_cov_demi(np.eye(control.x.size) * 0.5)
    control.x[0:5] = [2.0, 0.0, -2.0, -2.0, 1.0]

    # Compute cost and gradients
    cost, grad = dassflow1d.calc_cost_and_gradients(model, control, obs)
    cost_target = 27.451394309121454
    #grad_target = 0.5 * np.array([5.17981200e+01, 4.14673071e-01, 7.00174593e-04, 5.71318712e-01, 5.73441795e+00])
    grad_target = 0.5 * np.array([5.17981248e+01, 4.14673721e-01, 7.00159669e-04, 5.71388443e-01, 5.73614003e+00])
    print('cost',cost,cost_target)
    print('grad cast', grad, grad_target)
    assert(np.allclose(cost, cost_target))
    assert(np.allclose(grad, grad_target))


#def gradient_test(dt):

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
    #model.dt = dt
    #model.theta_preissmann = 0.9
    #model.set_scheme("preissmann")

    ## Create observations
    #obs = m_obs.all_observed([0.0, 1800.0, 3600.0, 5400.0, 7200.0], mesh)
    
    ## Time loop
    #model.generate_observations(mesh, obs)
    
    ## Setup control
    #control = m_control.Control()
    #control.add_bc_in_control(model, 0)
    #control.x[0:5] = [2.0, 2.0, 4.0, 2.0, 1.5]
    #diff_x = np.ones(5) * 0.1
    
    ## Compute initial cost and gradients
    ##cost0 = dassflow1d.calc_cost(model, mesh, imp, dof, control, obs)
    
    #cost0, grad0 = dassflow1d.calc_cost_and_gradients(model, mesh, imp, dof, control, obs)
    #print("cost0=", cost0)
    #print("grad0=", grad0)
    
    #cost_diff = np.sum(diff_x * grad0)

    #grad = np.zeros(control.x.size)
    #x0 = control.x.copy()
    #for eps in [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001]:
        #control.x[:] = x0 + eps * diff_x
        #cost = dassflow1d.calc_cost(model, mesh, imp, dof, control, obs)
        ##control.x[i] = x0
        ##print("(%i):cost=" % i, cost)
        ##grad[i] = (cost - cost0) / 0.001
        #print("eps=%.8f:%.8e" % (eps, 1.0 - (cost - cost0) / (eps * cost_diff)))
    

# To run tests without pytest (debug)
if __name__ == "__main__":
    test_all_observed()
    test_all_observed_change_of_variable()
    #gradient_test(dt=sys.argv[1])
