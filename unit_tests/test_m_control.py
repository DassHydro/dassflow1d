import numpy as np
import os

import dassflow1d
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_control as m_control


def test_control_setup():

    # Test initialisation of a unknowns on a mesh with two xs
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    model = m_sw_mono.Model(mesh)
    
    model.bc[0].id = "discharge"
    model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    model.bc[1].id = "elevation"
    model.bc[1].ts = m_sw_mono.Timeseries(0, 0.1)
    
    control = m_control.Control()
    control.add_bc_in_control(model, 0)
    
    print(control.items)
    for item in control.items:
        print(item)
    print(control.x)
    
    assert(control.items[0].id == b"BC001 ")
    assert(control.items[0].nx == 5)
    assert(np.allclose(control.x, [1.0, 2.0, 5.0, 3.0, 1.0]))
    #assert(unk.q.size == 2)
    #assert(unk.h.size == 2)
    #assert(unk.sg.size == 2)
    #assert(unk.sf.size == 2)


def test_apply_control():

    # Test initialisation of a unknowns on a mesh with two xs
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    model = m_sw_mono.Model(mesh)
    
    model.bc[0].id = "discharge"
    model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    model.bc[1].id = "elevation"
    model.bc[1].ts = m_sw_mono.Timeseries(0, 0.1)
    
    control = m_control.Control()
    control.add_bc_in_control(model, 0)
    
    control.x[0:5] = 1.0 
    
    control.apply_control(model)
    
    print(control.items[0])
    print(control.x)
    
    print(model.bc[0].ts)
    assert(control.items[0].id == b"BC001 ")
    assert(control.items[0].nx == 5)
    assert(np.allclose(control.x, 1.0))
    
    #assert(unk.a.size == 2)
    #assert(unk.q.size == 2)
    #assert(unk.h.size == 2)
    #assert(unk.sg.size == 2)
    #assert(unk.sf.size == 2)


def test_change_of_variable():

    # Test initialisation of a unknowns on a mesh with two xs
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    model = m_sw_mono.Model(mesh)
    
    model.bc[0].id = "discharge"
    model.bc[0].set_timeseries(t=[0.0, 1800.0, 3600.0, 5400.0, 7200.0], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    model.bc[1].id = "elevation"
    model.bc[1].ts = m_sw_mono.Timeseries(0, 0.1)
    
    control = m_control.Control()
    control.add_bc_in_control(model, 0)
    control.set_prior_cov_demi(np.eye(control.x.size) * 0.5)
    control.x[:] = 0.1
    
    control.apply_control(model)
    
    print(control.items[0])
    print(control.x)
    print(control.bdemi)
    
    print(model.bc[0].ts)
    assert(control.items[0].id == b"BC001 ")
    assert(control.items[0].nx == 5)
    assert(np.allclose(control.x, 0.1))
    assert(np.allclose(model.bc[0].ts.y, [1.05, 2.05, 5.05, 3.05, 1.05]))
    
    #assert(unk.a.size == 2)
    #assert(unk.q.size == 2)
    #assert(unk.h.size == 2)
    #assert(unk.sg.size == 2)
    #assert(unk.sf.size == 2)


# To run tests without pytest (debug)
if __name__ == "__main__":
    #test_control_setup()
    #test_apply_control()
    test_change_of_variable()
