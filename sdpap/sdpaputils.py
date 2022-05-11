#!/usr/bin/env python
"""sdpaputils.py

This is the module of sdpap

December 2010, Kenta KATO
"""

__all__ = ['get_dualitygap', 'get_primalerror', 'get_dualerror']

from .symcone import SymCone
from scipy.sparse.linalg.eigen import arpack
from scipy.sparse import csc_matrix
from scipy import sparse

def get_dualitygap(x, y, b, c):
    """Get duality gap of result

    Args:
      x: A primal result
      y: A dual result
      b: A constraint vector of primal
      c: A constraint vector of dual

    Returns:
      A float value of duality gap
    """

    if not sparse.isspmatrix_csc(x):
        x = x.tocsc()
    if not sparse.isspmatrix_csc(y):
        y = y.tocsc()
    if not sparse.isspmatrix_csc(b):
        b = b.tocsc()
    if not sparse.isspmatrix_csc(c):
        c = c.tocsc()

    objP = (c.T * x)[0,0]
    objD = (b.T * y)[0,0]

    return abs(objP - objD) / max(1.0, (abs(objP) + abs(objD)) / 2)


def get_primalerror(x, A, b, J):
    """Get primal feasible error

    Args:
      x: A primal result
      A: A constraint matrix
      b: A constraint vector of primal
      J: A SymCone

    Returns:
      A float value of primal feasible error
    """
    if not sparse.isspmatrix_csc(x):
        x = x.tocsc()
    if not sparse.isspmatrix_csc(A):
        A = A.tocsc()
    if not sparse.isspmatrix_csc(b):
        b = b.tocsc()

    delta = A * x - b
    maxerr = 0.0

    # Feasible error for K.f
    if J.f > 0:
        err = max(abs(delta[0:J.f, :]))[0, 0]
        maxerr = max(maxerr, err)

    if J.l > 0:
        err = min(delta[J.f:(J.f + J.l), :])[0,0]
        if err < 0:
            maxerr = max(maxerr, abs(err))

    if len(J.q) > 0:
        offset = J.f + J.l
        for k in J.q:
            vec = delta[offset:(offset + k), :]
            vec0 = vec[0, 0]
            vec1 = vec[1:, 0]
            err = vec0 ** 2 - (vec1.T * vec1)[0,0]
            if err < 0:
                maxerr = max(maxerr, abs(err))

            offset += k

    if len(J.s) > 0:
        offset = J.f + J.l + sum(J.q)
        s_err = 0.0
        for k in J.s:
            vec = delta[offset:(offset + k ** 2), :]
            mat_row = [i // k for i in vec.indices]
            mat_col = [i % k for i in vec.indices]
            mat_val = vec.data
            mat = csc_matrix((mat_val, (mat_row, mat_col)), shape=(k, k))
            # scipy.sparse.linalg.eigs raises a TypeError if mat is sparse and
            # k (1 in this case) is >= n (no. of rows of mat) - 1
            if mat.shape[0] <= 2:
                eig = arpack.eigs(mat.toarray(), k=1, which='SM',
                                  return_eigenvectors=False)[0]
            else:
                eig = arpack.eigs(mat, k=1, which='SM',
                                  return_eigenvectors=False)[0]

            if eig < 0:
                maxerr = max(maxerr, abs(eig))

            offset += k ** 2

    return maxerr


def get_dualerror(y, A, c, K):
    """Get dual feasible error

    Args:
      y: A dual result
      A: A constraint matrix
      c: A constraint vector of dual
      K: A SymCone

    Returns:
      A float value of primal feasible error
    """
    return get_primalerror(y, -A.T, -c, K)
