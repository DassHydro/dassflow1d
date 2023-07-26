import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.optimize

import dassflow1d
import dassflow1d.m_control as m_control
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_obs as m_obs


def test_forward_case_discharge():

    # Read mesh
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    mesh.set_strickler_type("powerlaw_h")
    mesh.set_uniform_strickler_parameters([20.0, 0.0])
    #mesh.set_bathy_field(x=)    
    
    # Create model
    model = m_sw_mono.Model(mesh)
    model.bc[0].id = "discharge"
    model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    model.bc[1].id = "normal_depth"
    model.ts = 0
    model.te = 7200
    model.dt = 60
    model.dtout = -1
    model.theta_preissmann = 0.667
    #model.set_scheme("implicit_diffusive_wave")
    model.set_scheme("preissmann")

    # Create observations
    obs = m_obs.all_observed([0.0, 1800.0, 3600.0, 5400.0, 7200.0], mesh)
    
    # Generate observations
    model.generate_observations(obs)
    
    # Setup control for optimisation
    control = m_control.Control()
    control.add_bc_in_control(model, 0)
    x_target = control.x.copy()
    control.x[0:5] = [2.0, 2.0, 6.0, 2.0, 1.5]

    # Perform optimisation
    res = scipy.optimize.minimize(dassflow1d.calc_cost_and_gradients_scipy_minimize, control.x, 
                                  args=(model, control, obs), jac=True,
                                  method='L-BFGS-B', options={"disp": False})
    print('res',res.x)
    print('target',x_target)
    print('jac',res.jac)
    print(obs.obs)
    assert(np.allclose(res.x, x_target, atol=1e-6))


#####def test_forward_case_bathy(verbose=False, plot=False):

    ###### Read mesh
    #####mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    #####mesh.set_uniform_strickler_parameters([20.0, 0.0])
    ######mesh.set_bathy_field(x=)    
    
    ###### Create model and unknowns
    #####model = m_model.Model(mesh)
    #####dof = m_model.Unknowns(mesh)
    
    ###### Setup boundary conditions
    #####model.bc[0].id = "discharge"
    ######model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    #####model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[0.5, 1.0, 2.5, 1.5, 0.5])
    #####model.bc[1].id = "normal_depth"
    #####imp = m_model.Implicitmatrix(model, mesh)
    #####model.ts = 0
    #####model.te = 7200
    #####model.dt = 60
    #####model.dtout = -1
    #####model.theta_preissmann = 0.667

    ###### Generate observations
    #####obs = m_obs.all_observed([0.0, 1800.0, 3600.0, 5400.0, 7200.0], mesh)
    #####model.generate_observations(mesh, imp, dof, obs)
    
    ###### Setup control for optimisation
    #####control = m_control.Control()
    #####control.add_bathy_in_control(model, mesh)
    #####x_target = control.x.copy()
    #####control.x += np.random.normal(0.0, 0.05, control.x.size)
    ######control.x[:] = 0.1
    #####control.x[:] = [1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    #####x_prior = control.x.copy()
    
    #####cost, grad = dassflow1d.calc_cost_and_gradients_scipy_minimize(x_prior, model, mesh, imp, dof, control, obs)
    ######print ("cost=", cost)
    ######print ("grad=", grad)
    #####if plot:
      #####plt.plot(obs.obs[0, :], obs.est[0, :], 'k.')
      #####plt.show()
    
    ######cost = dassflow1d.calc_cost(model, mesh, imp, dof, control, obs)
    ######if plot:
      ######for station in obs.stations:
        ######plt.plot(obs.obs[0, station.offset:station.offset+station.t.size], 'r.')
        ######plt.plot(obs.est[0, station.offset:station.offset+station.t.size], 'b-')
        ######plt.show()
      ######plt.plot(obs.obs[0, :], obs.est[0, :], 'k.')
      ######plt.show()
    
    ###### Perform optimisation
    #####res = scipy.optimize.minimize(dassflow1d.calc_cost_and_gradients_scipy_minimize, control.x, 
                                  #####args=(model, mesh, imp, dof, control, obs), jac=True,
                                  #####method='L-BFGS-B', options={"disp": verbose},
                                  #####callback=callback_minimize)
    #####if res.fun > 0.1:
      #####return
    #####if plot:
      
      #####x_station = []
      #####for station in obs.stations:
        #####x = 0
        #####for ics in station.ics:
          #####x += mesh.cs[ics-1].x
        #####x_station.append(x / station.ics.size)
      #####x_station = np.array(x_station)
      
      ######plt.plot(x_station, x_target, 'r-')
      ######plt.plot(x_station, x_prior, 'k--')
      ######plt.plot(x_station, res.x, 'b-')
      ######plt.show()

      #####H_obs = np.zeros((x_station.size, 5))
      #####H_est = np.zeros((x_station.size, obs.stations[0].t.size))
      #####ista = 0
      #####for station in obs.stations:
        #####H_obs[ista, :] =  obs.obs[0, station.offset:station.offset+station.t.size]
        #####H_est[ista, :] =  obs.est[0, station.offset:station.offset+station.t.size]
        #####ista += 1
        ######plt.plot(obs.obs[0, station.offset:station.offset+station.t.size], 'r.')
        ######plt.plot(obs.est[0, station.offset:station.offset+station.t.size], 'b-')
        ######plt.show()
      #####plt.plot(x_station, x_target, 'r-')
      #####plt.plot(x_station, x_prior, 'k--')
      #####plt.plot(x_station, res.x, 'g-')
      #####plt.plot(x_station, H_obs[:, 0], 'r.')
      #####plt.plot(x_station, H_est[:, 0], 'b-')
      #####plt.show()
        
      #####plt.plot(obs.obs[0, :], obs.est[0, :], 'k.')
      #####plt.show()
      
    ######print(obs.est)

    ######assert(np.allclose(res.x, x_target, atol=1e-6))


# To run tests without pytest (debug)
if __name__ == "__main__":
    test_forward_case_discharge()
