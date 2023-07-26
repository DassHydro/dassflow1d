import numpy as np
import os

import dassflow1d.m_mesh as m_mesh
import dassflow1d


def test_read_spatial_field_01():

    # Test reading of mesh01.geo
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh01.geo"))

    # Test reading of alpha and beta defined at every cross-sections
    mesh.read_spatial_field(os.path.join("data", "strickler_mesh01.dat"))
    
    # Test values
    css = [mesh.cs[i] for i in range(mesh.seg[0].first_cs-1, mesh.seg[0].last_cs)]
    for i in range(0, 11):
        assert(np.isclose(css[i].strickler_params[0], (i+1)*0.01))
        assert(np.isclose(css[i].strickler_params[1], (i+1)*0.001))


if __name__ == "__main__":
  test_read_spatial_field_01()
