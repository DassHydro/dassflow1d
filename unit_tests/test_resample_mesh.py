import numpy as np
import os

import dassflow1d.m_mesh as m_mesh
import dassflow1d


#def test_resample_mesh_01():

    ## Creation du maillage de base (2 cross-sections, 1 segment)
    #mesh = m_mesh.Mesh(6, 1)
    
    #h = np.array([1.0, 2.0, 3.0, 4.0])
    #w = np.array([100.0, 200.0, 300.0, 400.0])
    #b = np.array([0.5, 0.4, 0.3, 0.2, 0.1, 0.0])
    
    #for i in range(0, 6):
        #mesh.cs[i].set_levels(h+b
        
    #xs = m_mesh.Crosssection(2)
    #assert(xs.level_heights.size == 2)
    #assert(xs.level_widths.size == 2)
    #assert(xs.poly.shape == (2, 2))
    #assert(np.allclose(xs.level_heights, [0.0, 0.0]))
    #assert(np.allclose(xs.level_widths, [0.0, 0.0]))

    ## Test initialisation of a cross-section with shape_model 'linear'
    #xs = m_mesh.Crosssection(2, 'linear')
    #assert(xs.level_heights.size == 2)
    #assert(xs.level_widths.size == 2)
    #assert(xs.poly.shape == (2, 2))
    #assert(np.allclose(xs.level_heights, [0.0, 0.0]))
    #assert(np.allclose(xs.level_widths, [0.0, 0.0]))

    ## Test initialisation of a cross-section with shape_model 'cubic_spline'
    #xs = m_mesh.Crosssection(5, 'cubic_spline')
    #assert(xs.level_heights.size == 5)
    #assert(xs.level_widths.size == 5)
    #assert(xs.poly.shape == (4, 5))
    #assert(np.allclose(xs.level_heights, [0.0, 0.0, 0.0, 0.0, 0.0]))
    #assert(np.allclose(xs.level_widths, [0.0, 0.0, 0.0, 0.0, 0.0]))


def test_Mesh_initialise():

    # Test initialisation of a mesh with two xs
    mesh = m_mesh.Mesh(2)
    assert(mesh.ncs == 2)


def test_Mesh_read_01():

    # Test reading of mesh01.geo
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh01.geo"))
    
    # Assert number of cross-sections and segments
    assert(mesh.ncs == 17)
    assert(mesh.nseg == 1)
    
    # Assert cross-sections
    x = np.zeros(17)
    x[0:2] = 0.0
    x[2:-2] = np.linspace(0.0,  300.0, 13, endpoint=True)
    x[-2:] = 300.0
    bathy_exact = (300-x) * 0.01
    area_cum_exact = 20000.0 - (300-x) * 0.01
    perim_cum_exact = (20000.0 - (300-x) * 0.01) * 2.0 + 1.0
    for ics in range(0, 17):
      cs = mesh.cs[ics]
      assert(cs.nlevels == 1)
      assert(cs.bathy == bathy_exact[ics])
      assert(cs.level_heights[:] == [20000.0])
      assert(cs.level_widths[:] == [1.0])
      assert(np.allclose(cs.poly, [[0.0], [1.0]]))
      assert(np.allclose(cs.area_cum[1], [area_cum_exact[ics]]))
      assert(np.allclose(cs.perim_cum[1], [perim_cum_exact[ics]]))
    
    # Assert segment
    assert(mesh.seg[0].first_cs == 3)
    assert(mesh.seg[0].last_cs == 15)
    assert(mesh.seg[0].us_seg == [-1])
    assert(mesh.seg[0].ds_seg == -1)


def test_Mesh_read_02():

    # Test reading of mesh01.geo
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh02.geo"))
    print(mesh)
    
    # Assert number of cross-sections and segments
    assert(mesh.ncs == 30)
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
    area_cum_exact = (1000.0 - bathy_exact) * 10.0
    perim_cum_exact = (1000.0 - bathy_exact) * 2.0 + 10.0
    for ics in range(0, 10):
      cs = mesh.cs[ics]
      assert(cs.nlevels == 1)
      assert(cs.bathy == bathy_exact[ics])
      assert(cs.level_heights[:] == [1000.0])
      assert(cs.level_widths[:] == [10.0])
      assert(np.allclose(cs.poly, [[0.0], [10.0]]))
      assert(np.allclose(cs.area_cum[1], [area_cum_exact[ics]]))
      assert(np.allclose(cs.perim_cum[1], [perim_cum_exact[ics]]))
      #assert(np.allclose(cs.area_cum, [area_cum_exact[ics]]))
      #assert(np.allclose(cs.perim_cum, [perim_cum_exact[ics]]))
    
    # Assert cross-sections for segment 2
    x = np.zeros(10)
    x[0:2] = 0.0
    x[2:-2] = np.linspace(0.0,  500.0, 6, endpoint=True)
    x[-2:] = 500.0
    bathy_exact = 0.5 + (500-x) * 0.001
    area_cum_exact = (1000.0 - bathy_exact) * 10.0
    perim_cum_exact = (1000.0 - bathy_exact) * 2.0 + 10.0
    for ics in range(10, 20):
      cs = mesh.cs[ics]
      assert(cs.nlevels == 1)
      assert(cs.bathy == bathy_exact[ics-10])
      assert(cs.level_heights[:] == [1000.0])
      assert(cs.level_widths[:] == [10.0])
      assert(np.allclose(cs.poly, [[0.0], [10.0]]))
      assert(np.allclose(cs.area_cum[1], [area_cum_exact[ics-10]]))
      assert(np.allclose(cs.perim_cum[1], [perim_cum_exact[ics-10]]))
      #assert(np.allclose(cs.area_cum, [area_cum_exact[ics-10]]))
      #assert(np.allclose(cs.perim_cum, [perim_cum_exact[ics-10]]))
    
    # Assert cross-sections for segment 2
    x = np.zeros(10)
    x[0:2] = 0.0
    x[2:-2] = np.linspace(0.0,  500.0, 6, endpoint=True)
    x[-2:] = 500.0
    bathy_exact = (500-x) * 0.001
    area_cum_exact = (1000.0 - bathy_exact) * 10.0
    perim_cum_exact = (1000.0 - bathy_exact) * 2.0 + 10.0
    for ics in range(20, 30):
      cs = mesh.cs[ics]
      assert(cs.nlevels == 1)
      assert(cs.bathy == bathy_exact[ics-20])
      assert(cs.level_heights[:] == [1000.0])
      assert(cs.level_widths[:] == [10.0])
      assert(np.allclose(cs.poly, [[0.0], [10.0]]))
      assert(np.allclose(cs.area_cum[1], [area_cum_exact[ics-20]]))
      assert(np.allclose(cs.perim_cum[1], [perim_cum_exact[ics-20]]))
      #assert(np.allclose(cs.area_cum, [area_cum_exact[ics-20]]))
      #assert(np.allclose(cs.perim_cum, [perim_cum_exact[ics-20]]))



def test_Mesh_read_03():

    # Test reading of mesh01.geo
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh03.geo"))
    
    # Assert number of cross-sections and segments
    assert(mesh.ncs == 30)
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
    for ics in range(0, 10):
      cs = mesh.cs[ics]
      assert(cs.nlevels == 3)
      assert(cs.bathy == bathy_exact[ics])
      assert(np.allclose(cs.level_heights, [bathy_exact[ics]+1.0, bathy_exact[ics]+2.0, bathy_exact[ics]+3.0]))
      assert(np.allclose(cs.level_widths, [100.0, 150.0, 1000.0]))
      assert(cs.ob_levels[0] == 3)
      assert(cs.ob_levels[1] == 2)
      poly_b_exact = cs.level_widths[:] - poly_a_exact[:] * cs.level_heights[:]
      assert(np.allclose(cs.poly, [poly_a_exact, poly_b_exact]))
      assert(np.allclose(cs.area_cum, area_cum_exact))
      assert(np.allclose(cs.perim_cum, perim_cum_exact))
    
    # Assert cross-sections for segment 2
    x = np.zeros(10)
    x[0:2] = 0.0
    x[2:-2] = np.linspace(0.0,  500.0, 6, endpoint=True)
    x[-2:] = 500.0
    bathy_exact = 0.5 + (500-x) * 0.001
    area_cum_exact = [0.0, 100.0, 225.0, 800.0]
    dperim = [102.0, 2*np.sqrt(25.0**2+1), 2*np.sqrt(425.0**2+1)]
    perim_cum_exact = [100.0, dperim[0], np.sum(dperim[0:2]), np.sum(dperim)]
    poly_a_exact = [0.0, 50.0, 850.0]
    for ics in range(10, 20):
      cs = mesh.cs[ics]
      assert(cs.nlevels == 3)
      assert(cs.bathy == bathy_exact[ics-10])
      assert(np.allclose(cs.level_heights, [bathy_exact[ics-10]+1.0, bathy_exact[ics-10]+2.0, bathy_exact[ics-10]+3.0]))
      assert(np.allclose(cs.level_widths, [100.0, 150.0, 1000.0]))
      assert(cs.ob_levels[0] == 3)
      assert(cs.ob_levels[1] == 2)
      poly_b_exact = cs.level_widths[:] - poly_a_exact[:] * cs.level_heights[:]
      assert(np.allclose(cs.poly, [poly_a_exact, poly_b_exact]))
      assert(np.allclose(cs.area_cum, area_cum_exact))
      assert(np.allclose(cs.perim_cum, perim_cum_exact))
    
    # Assert cross-sections for segment 2
    x = np.zeros(10)
    x[0:2] = 0.0
    x[2:-2] = np.linspace(0.0,  500.0, 6, endpoint=True)
    x[-2:] = 500.0
    bathy_exact = (500-x) * 0.001
    area_cum_exact = [0.0, 100.0, 225.0, 800.0]
    dperim = [102.0, 2*np.sqrt(25.0**2+1), 2*np.sqrt(425.0**2+1)]
    perim_cum_exact = [100.0, dperim[0], np.sum(dperim[0:2]), np.sum(dperim)]
    poly_a_exact = [0.0, 50.0, 850.0]
    for ics in range(20, 30):
      cs = mesh.cs[ics]
      assert(cs.nlevels == 3)
      assert(cs.bathy == bathy_exact[ics-20])
      assert(np.allclose(cs.level_heights, [bathy_exact[ics-20]+1.0, bathy_exact[ics-20]+2.0, bathy_exact[ics-20]+3.0]))
      assert(np.allclose(cs.level_widths, [100.0, 150.0, 1000.0]))
      assert(cs.ob_levels[0] == 3)
      assert(cs.ob_levels[1] == 2)
      poly_b_exact = cs.level_widths[:] - poly_a_exact[:] * cs.level_heights[:]
      assert(np.allclose(cs.poly, [poly_a_exact, poly_b_exact]))
      assert(np.allclose(cs.area_cum, area_cum_exact))
      assert(np.allclose(cs.perim_cum, perim_cum_exact))



def test_Mesh_apply_bathy_field_01():

    # Read of mesh01.geo
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh01.geo"))
    
    # Setup bathy field
    mesh.set_bathy_field_linear([0.0, 300.0], [0.0, 1.0])
    mesh.apply_bathy_field(0.1)
    bathy = np.zeros(13)
    for ics in range(2, 15):
      bathy[ics-2] = mesh.cs[ics].bathy


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
      ics += 1


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


if __name__ == "__main__":
  #test_Crosssection_initialise()
  #test_Mesh_initialise()
  #test_Mesh_read_01()
  #test_Mesh_read_02()
  test_Mesh_read_03()
  #test_Mesh_apply_bathy_field_01()
  #test_strickler_fields_k_constant()
  #test_strickler_fields_k_powerlaw_h()
