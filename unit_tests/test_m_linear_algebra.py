import numpy as np
import os

import dassflow1d
import dassflow1d.m_linear_algebra as m_linear_algebra


def test_matrix_csr():

    # Test initialisation of a unknowns on a mesh with two xs
    array = np.array([[1.0, 0.0, 2.0], [0.0, 3.0, 4.0], [5.0, 6.0, 0.0]])
    M = m_linear_algebra.matrixcsr_from_numpy_array(array)

    assert(M.n == 3)
    assert(M.nnz == 6)
    assert(np.all(M.irow == [0, 0, 1, 1, 2, 2]))
    assert(np.all(M.icol == [0, 2, 1, 2, 0, 1]))
    assert(np.allclose(M.anz, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))


# To run tests without pytest (debug)
if __name__ == "__main__":
    test_matrix_csr()
    #test_apply_control()
