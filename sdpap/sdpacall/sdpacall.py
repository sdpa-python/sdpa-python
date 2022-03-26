#!/usr/bin/env python
"""sdpacall.py

A module of sdpap
Convert to MatData and call sdpa

December 2010, Kenta KATO
"""

__all__ = ['solve_sdpa']

from . import sdpa
from sdpap.matdata import MatData
from sdpap.symcone import SymCone
from scipy.sparse import csr_matrix
from scipy import sparse, array


def solve_sdpa(A, b, c, K, option):
    """Solve SDP with sdpa

    Args:
      A, b, c: Sparse matrix which denotes SDP problem
      K: SymCone
      option: Dictionary of option parameters

    Returns:
      A tuple (x, y, info). Each members are:
        x: Primal result
        y: Dual result
        info: Result SDPA information
    """
    At = (A.T).tocsc()
    At.sort_indices()
    if not sparse.isspmatrix_csc(b):
        b = b.tocsc()
    if not sparse.isspmatrix_csc(c):
        c = c.tocsc()
    b.sort_indices()
    c.sort_indices()
    data_At = MatData(At)
    data_b = MatData(b)
    data_c = MatData(c)
    dict_K = K.todict()

    data_x, data_y, data_s, info = sdpa.sedumiwrap(data_At, data_b, data_c,
                                                   dict_K, option)

    x = csr_matrix(array(data_x)).T
    y = csr_matrix(array(data_y)).T
    s = csr_matrix(array(data_s)).T

    return x, y, s, info

