import numpy as np
#import matplotlib.pyplot as plt

import dassflow1d
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono


def test_Unknowns():

    # Test initialisation of a unknowns on a mesh with two xs
    mesh = m_mesh.Mesh(2)
    unk = m_sw_mono.Unknowns(mesh)
    assert(unk.a.size == 2)
    assert(unk.q.size == 2)
    assert(unk.h.size == 2)
    assert(unk.sg.size == 2)
    assert(unk.sf.size == 2)


# To run tests without pytest (debug)
if __name__ == "__main__":
    test_Crosssection_trapezoids()
