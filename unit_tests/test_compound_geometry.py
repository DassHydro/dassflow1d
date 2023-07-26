import numpy as np
#import matplotlib.pyplot as plt

import dassflow1d
import dassflow1d.m_mesh as m_mesh


def test_Crosssection_trapezoids():

    # Test initialisation of a cross-section with default shape_model
    cs = m_mesh.Crosssection(3)
    cs.bathy = 0.0
    cs.level_heights[:] = [1.0, 2.0, 3.0]
    cs.level_widths[:] = [10.0, 100.0, 120.0]
    cs.ob_levels[:] = [1, 2]
    cs.update_geometry()
    assert(cs.level_heights.size == 3)
    assert(cs.level_widths.size == 3)
    assert(cs.poly.shape == (2, 3))
    assert(np.allclose(cs.level_heights, [1.0, 2.0, 3.0]))
    assert(np.allclose(cs.level_widths, [10.0, 100.0, 120.0]))
    
    Z = np.linspace(0, 3.0, 7, endpoint=True)
    W = np.zeros(Z.size)
    Wexact = np.array([10.0, 10.0, 10.0, 55.0, 100.0, 110.0, 120.0])
    for i, h in enumerate(Z):
        cs.update_level(h)
        height = cs.bathy + h
        width = cs.width(h)
        W[i] = width
        
    assert(np.allclose(W, Wexact))

    A = np.zeros((Z.size, 3))
    areas = np.zeros(3)
    Aexact = np.array([[0.0, 0.0, 0.0], 
                       [0.0, 5.0, 0.0],
                       [0.0, 10.0, 0.0],
                       [5.625, 20.625, 0.0],
                       [22.5, 42.5, 0.0],
                       [46.25, 70.0, 1.25],
                       [72.5, 97.5, 5.0]])
    for i, h in enumerate(Z):
        cs.update_level(h)
        height = cs.bathy + h
        cs.areas_compound_channel(h, areas)
        A[i, :] = areas[:]
    print(A)
        
    assert(np.allclose(A, Aexact))

    P = np.zeros((Z.size, 3))
    perims = np.zeros(3)
    Pexact = np.array([[0.0, 10.0, 0.0],
                       [0.0, 11.0, 0.0],
                       [0.0, 12.0, 0.0],
                       [22.50555487, 34.50555487, 0.0],
                       [45.01110974, 57.01110974, 0.0],
                       [50.03604755, 57.01110974, 5.02493781],
                       [55.06098536, 57.01110974, 10.04987562]])
    for i, h in enumerate(Z):
        cs.update_level(h)
        height = cs.bathy + h
        cs.perimeters_compound_channel(h, perims)
        P[i, :] = perims
    print(P)
        
    assert(np.allclose(P, Pexact))


def test_Crosssection_quadratic():
  
    # TODO
    pass


def test_Crosssection_cubic():
  
    # TODO
    pass


# To run tests without pytest (debug)
if __name__ == "__main__":
    test_Crosssection_trapezoids()
