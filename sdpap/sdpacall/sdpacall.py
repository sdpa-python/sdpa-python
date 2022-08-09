#!/usr/bin/env python
"""
Routines to convert to MatData and call SDPA (C API)
This file is a component of SDPAP
Copyright (C) 2010-2022 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

December 2010: Originally written by Kenta Kato
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

