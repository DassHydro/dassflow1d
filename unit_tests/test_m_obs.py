import numpy as np
import os

import dassflow1d
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_obs as m_obs


def test_Obsstation():

    # Test initialisation of a unknowns on a mesh with two xs
    obsstation = m_obs.Obsstation([1,2,3], [0.0, 60.0, 120.0], [[1.0, 1.0, 1.0],[0.1, 0.1, 0.1]])
    assert(np.allclose(obsstation.t, [0.0, 60.0, 120.0]))
    assert(np.allclose(obsstation.obs, [[1.0, 1.0, 1.0],[0.1, 0.1, 0.1]]))


def test_all_observed():

    # Test initialisation of a unknowns on a mesh with two xs
    mesh = dassflow1d.read_mesh(os.path.join("data", "mesh_case_forward_channel.geo"))
    obs = m_obs.all_observed([0.0, 60.0, 120.0], mesh)
    #print(obs)
    #print(obs.stations[0])
    ista = 0
    for station in obs.stations:
      assert(np.allclose(station.t, [0.0, 60.0, 120.0]))
      assert(station.offset == ista * 3)
      ista += 1


# To run tests without pytest (debug)
if __name__ == "__main__":
    test_Obsstation()
    test_all_observed()
