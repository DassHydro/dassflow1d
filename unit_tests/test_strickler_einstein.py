import numpy as np
import matplotlib.pyplot as plt
import os

import dassflow1d.m_mesh as m_mesh
import dassflow1d


def test_strickler_einstein_mesh03():

    # Load mesh mesh03.geo
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh03.geo"))

    # Set strickler type
    mesh.set_strickler_type("Einstein")
    mesh.set_uniform_strickler_parameters([10.0, 20.0, 15.0])
    
    # Assert strickler type code (mesh side)
    assert(mesh.strickler_type_code == m_mesh.strickler_type_einstein)
    
    # Test strickler values
    assert(mesh.nseg == 3)
    
    # Assert segments
    assert(mesh.seg[0].first_cs == 3)
    assert(mesh.seg[0].last_cs == 8)
    assert(mesh.seg[0].us_seg == [-1])
    assert(mesh.seg[0].ds_seg == 2)
    assert(mesh.seg[1].first_cs == 13)
    assert(mesh.seg[1].last_cs == 18)
    assert(mesh.seg[1].us_seg == [1])
    assert(mesh.seg[1].ds_seg == 3)
    assert(mesh.seg[2].first_cs == 23)
    assert(mesh.seg[2].last_cs == 28)
    assert(mesh.seg[2].us_seg == [2])
    assert(mesh.seg[2].ds_seg == -1)
    
    # Assert cross-sections for segment 1
    x = np.zeros(10)
    x[0:2] = 0.0
    x[2:-2] = np.linspace(0.0,  500.0, 6, endpoint=True)
    x[-2:] = 500.0
    bathy_exact = 1.0 + (500-x) * 0.001
    area_cum_exact = [0.0, 100.0, 225.0, 800.0]
    dperim = [102.0, 2*np.sqrt(25.0**2+1), 2*np.sqrt(425.0**2+1)]
    perim_cum_exact = [100.0, dperim[0], np.sum(dperim[0:2]), np.sum(dperim)]
    poly_a_exact = [0.0, 50.0, 850.0]
    for ics in range(0, 1):
        
        # Assert strickler type code (cross-section side)
        cs = mesh.cs[ics]
        assert(cs.strickler_type_code == m_mesh.strickler_type_einstein)
        
        h = np.linspace(0.0, 4.0, 91, endpoint=True)
        Z = h + cs.bathy
        
        A = np.zeros(h.size)
        P = np.zeros(h.size)
        Q = np.zeros(h.size)
        
        for i in range(0, h.size):
            cs.update_level(h[i])
            K = cs.strickler(h[i])
            A[i] = cs.area(h[i])
            P[i] = cs.perimeter(h[i])
            U = K * (A[i]/P[i])**2./3. * np.sqrt(0.001)
            Q[i] = K * A[i] * (A[i]/P[i])**2./3. * np.sqrt(0.001)
            #print(i, h[i], K, A[i], P[i], A[i]/P[i], U, Q)
        
        plt.plot(-0.5 * cs.level_widths, cs.level_heights, 'b-')
        plt.plot([-0.5 * cs.level_widths[0]]*2, [cs.level_heights[0], cs.bathy], 'b-')
        plt.plot([-0.5 * cs.level_widths[0], 0.5 * cs.level_widths[0]], [cs.bathy, cs.bathy], 'b-')
        plt.plot([0.5 * cs.level_widths[0]]*2, [cs.level_heights[0], cs.bathy], 'b-')
        plt.plot(0.5 * cs.level_widths, cs.level_heights, 'b-')
        plt.plot(-0.5 * cs.level_widths[1], cs.level_heights[1], 'ro', label="banks")
        plt.plot(0.5 * cs.level_widths[1], cs.level_heights[1], 'ro')
        plt.xlabel("W (m)")
        plt.ylabel("H (m)")
        plt.legend()
        plt.show()

        plt.plot(h, Q, "b-")
        plt.xlabel("h (m)")
        plt.ylabel("Q (m3/s)")
        #plt.plot(h, A/P, "b-")
        plt.show()
        
      #assert(cs.nlevels == 3)
      #assert(cs.bathy == bathy_exact[ics])
      #assert(np.allclose(cs.level_heights, [bathy_exact[ics]+1.0, bathy_exact[ics]+2.0, bathy_exact[ics]+3.0]))
      #assert(np.allclose(cs.level_widths, [100.0, 150.0, 1000.0]))
      #assert(cs.ob_levels[0] == 3)
      #assert(cs.ob_levels[1] == 2)
      #poly_b_exact = cs.level_widths[:] - poly_a_exact[:] * cs.level_heights[:]
      #assert(np.allclose(cs.poly, [poly_a_exact, poly_b_exact]))
      #assert(np.allclose(cs.area_cum, area_cum_exact))
      #assert(np.allclose(cs.perim_cum, perim_cum_exact))
    
    ## Assert cross-sections for segment 2
    #x = np.zeros(10)
    #x[0:2] = 0.0
    #x[2:-2] = np.linspace(0.0,  500.0, 6, endpoint=True)
    #x[-2:] = 500.0
    #bathy_exact = 0.5 + (500-x) * 0.001
    #area_cum_exact = [0.0, 100.0, 225.0, 800.0]
    #dperim = [102.0, 2*np.sqrt(25.0**2+1), 2*np.sqrt(425.0**2+1)]
    #perim_cum_exact = [100.0, dperim[0], np.sum(dperim[0:2]), np.sum(dperim)]
    #poly_a_exact = [0.0, 50.0, 850.0]
    #for ics in range(10, 20):
      #cs = mesh.cs[ics]
      #assert(cs.nlevels == 3)
      #assert(cs.bathy == bathy_exact[ics-10])
      #assert(np.allclose(cs.level_heights, [bathy_exact[ics-10]+1.0, bathy_exact[ics-10]+2.0, bathy_exact[ics-10]+3.0]))
      #assert(np.allclose(cs.level_widths, [100.0, 150.0, 1000.0]))
      #assert(cs.ob_levels[0] == 3)
      #assert(cs.ob_levels[1] == 2)
      #poly_b_exact = cs.level_widths[:] - poly_a_exact[:] * cs.level_heights[:]
      #assert(np.allclose(cs.poly, [poly_a_exact, poly_b_exact]))
      #assert(np.allclose(cs.area_cum, area_cum_exact))
      #assert(np.allclose(cs.perim_cum, perim_cum_exact))
    
    ## Assert cross-sections for segment 2
    #x = np.zeros(10)
    #x[0:2] = 0.0
    #x[2:-2] = np.linspace(0.0,  500.0, 6, endpoint=True)
    #x[-2:] = 500.0
    #bathy_exact = (500-x) * 0.001
    #area_cum_exact = [0.0, 100.0, 225.0, 800.0]
    #dperim = [102.0, 2*np.sqrt(25.0**2+1), 2*np.sqrt(425.0**2+1)]
    #perim_cum_exact = [100.0, dperim[0], np.sum(dperim[0:2]), np.sum(dperim)]
    #poly_a_exact = [0.0, 50.0, 850.0]
    #for ics in range(20, 30):
      #cs = mesh.cs[ics]
      #assert(cs.nlevels == 3)
      #assert(cs.bathy == bathy_exact[ics-20])
      #assert(np.allclose(cs.level_heights, [bathy_exact[ics-20]+1.0, bathy_exact[ics-20]+2.0, bathy_exact[ics-20]+3.0]))
      #assert(np.allclose(cs.level_widths, [100.0, 150.0, 1000.0]))
      #assert(cs.ob_levels[0] == 3)
      #assert(cs.ob_levels[1] == 2)
      #poly_b_exact = cs.level_widths[:] - poly_a_exact[:] * cs.level_heights[:]
      #assert(np.allclose(cs.poly, [poly_a_exact, poly_b_exact]))
      #assert(np.allclose(cs.area_cum, area_cum_exact))
      #assert(np.allclose(cs.perim_cum, perim_cum_exact))
      


if __name__ == "__main__":
  test_strickler_einstein_mesh03()
