import numpy as np
#import matplotlib.pyplot as plt

import dassflow1d
import dassflow1d.m_mesh as m_mesh


def test_Crosssection_trapezoids():

    # Test initialisation of a cross-section with default shape_model
    cs = m_mesh.Crosssection(2)
    cs.bathy = 0.0
    cs.level_heights[:] = [1.0, 2.0]
    cs.level_widths[:] = [10.0, 100.0]
    cs.update_geometry()
    assert(cs.level_heights.size == 2)
    assert(cs.level_widths.size == 2)
    assert(cs.poly.shape == (2, 2))
    assert(np.allclose(cs.level_heights, [1.0, 2.0]))
    assert(np.allclose(cs.level_widths, [10.0, 100.0]))
    
    Z = np.linspace(0, 3.0, 7, endpoint=True)
    W = np.zeros(Z.size)
    Wexact = np.array([10.0, 10.0, 10.0, 55.0, 100.0, 145.0, 190.0])
    for i, h in enumerate(Z):
        cs.update_level(h)
        height = cs.bathy + h
        width = cs.width(h)
        W[i] = width
        
    assert(np.allclose(W, Wexact))

    A = np.zeros(Z.size)
    Aexact = np.array([0.0, 5.0, 10.0, 26.25, 65.0, 126.25, 210.0])
    for i, h in enumerate(Z):
        cs.update_level(h)
        height = cs.bathy + h
        area = cs.area(h)
        A[i] = area
        #print(i, h, area)
        
    assert(np.allclose(A, Aexact))

    P = np.zeros(Z.size)
    Pexact = np.zeros(Z.size)
    Pexact[0:3] = [10.0, 11.0, 12.0]
    Pexact[3:] = 12.0 + 2.0 * np.sqrt((Z[3:]-Z[2])**2 + ((Z[3:]-Z[2])*90.0/2.0)**2)
    for i, h in enumerate(Z):
        cs.update_level(h)
        height = cs.bathy + h
        perimeter = cs.perimeter(h)
        P[i] = perimeter
        
    assert(np.allclose(P, Pexact))

    h = np.zeros(Z.size)
    hexact = Z[:]
    for i, area in enumerate(A):
        h[i] = cs.atoh(area)
        print(area, h[i])
        
    assert(np.allclose(h, hexact))


def test_Crosssection_quadratic():
  
    # TODO
    pass


def test_Crosssection_cubic():
  
    # TODO
    pass


# To run tests without pytest (debug)
if __name__ == "__main__":
    test_Crosssection_trapezoids()
