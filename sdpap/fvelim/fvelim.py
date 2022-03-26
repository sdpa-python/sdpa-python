#!/usr/bin/env python
"""fvelim.py

A module of fvelim
Eliminate free variables

December 2010, Kenta KATO
"""

__all__ = ['split', 'eliminate', 'result_split', 'result_elimination']

from . import fvelimext
from sdpap.symcone import SymCone
from sdpap.matdata import MatData
from scipy.sparse import csc_matrix, csr_matrix, bmat, eye
from scipy.sparse.linalg import spsolve
from scipy import sparse, ones

def split(A, b, c, K, rho):
    """Eliminate free variables by Split method
    Args:
      A, b, c, K: SeDuMi format
      rho: Parameter of range in split method or
        pivoting in elimination method

    Returns:
      Converted problem (A2, b2, c2, K2).
    """
    # --------------------------------------------------
    # Type check
    # --------------------------------------------------
    if not sparse.isspmatrix_csc(A):
        A = A.tocsc()
    if not sparse.isspmatrix_csc(b):
        b = b.tocsc()
    if not sparse.isspmatrix_csc(c):
        c = c.tocsc()

    # --------------------------------------------------
    # Get K.f part
    # --------------------------------------------------
    A_f = A[:, 0:K.f]
    A_lqs = A[:, K.f:]
    c_f = c[0:K.f, :]
    c_lqs = c[K.f:, :]

    # ----------------------------------------
    # Split method
    # ----------------------------------------
    c_fp = c_f + csc_matrix(ones((c_f.shape[0], 1))) * rho
    c_fm =  -c_f + csc_matrix(ones((c_f.shape[0], 1))) * rho
    A_fp = A_f
    A_fm = -A_f
    A2 = bmat([[A_fp, A_fm, A_lqs]])
    b2 = b
    c2 = bmat([[c_fp],
               [c_fm],
               [c_lqs]])
    K2 = SymCone(f=0, l=2*K.f+K.l, q=K.q, s=K.s)

    return A2, b2, c2, K2


def eliminate(A, b, c, K, rho, zeroPoint):
    """Eliminate free variables by Elimination method

    Args:
      A, b, c, K: SeDuMi format
      rho: Parameter of range in split method or
        pivoting in elimination method
      zeroPoint: The zero point of matrix operation

    Returns:
      A Tuple (A2, b2, c2, K2, LiP, U, Q, LPA_B, LPb_B, cfQU, gamma).
      Each members are:
        A2, b2, c2, K2: Converted Problem
        LiP, U, Q, LPA_B, LPb_B, cfQU: For retrieving results.
          If method=='split', these value is None.
        gamma: Constant of objective value.
          If method=='split', gamma = 0.
        rank_Af: Matrix rank of K.f part of A
    """
    # --------------------------------------------------
    # Type check
    # --------------------------------------------------
    if not sparse.isspmatrix_csc(A):
        A = A.tocsc()
    if not sparse.isspmatrix_csc(b):
        b = b.tocsc()
    if not sparse.isspmatrix_csc(c):
        c = c.tocsc()

    # --------------------------------------------------
    # Get K.f part
    # --------------------------------------------------
    A_f = A[:, 0:K.f]
    A_lqs = A[:, K.f:]
    c_f = c[0:K.f, :]
    c_lqs = c[K.f:, :]

    # ----------------------------------------
    # Elimination method
    # ----------------------------------------
    if not sparse.isspmatrix_csc(A_f):
        A_f = A_f.tocsc()
    if not sparse.isspmatrix_csc(A_lqs):
        A_lqs = A_lqs.tocsc()

    A_f.sort_indices()
    A_lqs.sort_indices()
    data_Af = MatData(A_f)
    data_Alqs = MatData(A_lqs)
    data_LiP, data_U, data_V, Q, rank_Af = \
              fvelimext.lu(data_Af, data_Alqs, rho, zeroPoint)
    LiP = MatData(values=data_LiP[0],
                  rowind=data_LiP[1],
                  colptr=data_LiP[2],
                  size=(A_f.shape[0], A_f.shape[0])).tomatrix()
    U = MatData(values=data_U[0],
                rowind=data_U[1],
                colptr=data_U[2],
                size=(rank_Af, rank_Af)).tomatrix()
    if data_V != None:
        V = MatData(values=data_V[0],
                    rowind=data_V[1],
                    colptr=data_V[2],
                    size=(rank_Af, K.f - rank_Af)).tomatrix()

    # Make LPA, LPb
    LPA = LiP * A_lqs
    LPb = LiP * b

    # Devide to Basis and Nonbasis
    LPA_B = LPA[0:rank_Af, :]
    A2 = LPA[rank_Af:, :]
    LPb_B = LPb[0:rank_Af, :]
    b2 = LPb[rank_Af:, :]

    ordering = range(K.f)
    for i, index in enumerate(Q):
        tmp = ordering[i]
        ordering[i] = ordering[index]
        ordering[index] = tmp

    # Convert objective function
    cfQ = c_f[ordering]

    # Linear dependent case
    if rank_Af < K.f:
        cfQ1 = cfQ[0:rank_Af, :]
        cfQ2 = cfQ[rank_Af:, :]
        cfQ = cfQ1

    cfQU = csr_matrix(spsolve(U.T, cfQ))

    # Unbound check
    # Linear dependent case
    if rank_Af < K.f:
        cfQUV = cfQU * V
        submat = cfQ2 - cfQUV.T
        if submat.nnz > 0 and max(abs(submat.data)) > zeroPoint:
            print("%.15e" % max(abs(submat.data)))
            raise ValueError("This problem is unbound")

    cfQUA_B = cfQU * LPA_B
    c2 = c_lqs - cfQUA_B.T

    # Get gamma
    gamma = (cfQU * LPb_B).todense()[0,0]
    print("gamma = %f" % gamma)

    # Make K2
    K2 = SymCone(f=0, l=K.l, q=K.q, s=K.s)

    return A2, b2, c2, K2, LiP, U, Q, LPA_B, LPb_B, cfQU, gamma, rank_Af


def result_split(x2, y2, s2, K):
    """Get result of problem with free variables by Split method

    Args:
      x2, y2, s2: SDPA result
      K: SymCone of original problem

    Returns:
      Result of original problem (x, y, s)
    """
    x_fp = x2[0:K.f, :]
    x_fm = x2[K.f:(K.f + K.f), :]
    x_lqs = x2[(K.f + K.f):, :]
    x_f = x_fp - x_fm
    x = bmat([[x_f],
              [x_lqs]])

    s_lqs = s2[(K.f + K.f):, :]
    s = bmat([[csc_matrix((K.f, 1))],
              [s_lqs]])

    y = y2

    return x, y, s


def result_elimination(x2, y2, s2, K, LiP, U, Q, LPA_B, LPb_B, cfQU, rank_Af):
    """Get result of problem with free variables by Elimination method

    Args:
      x2, y2, s2: SDPA result
      K: SymCone of original problem
      LiP, U, Q, LPA_B, LPb_B, cfQU, rank_Af: Informations of LU factorization

    Returns:
      Result of original problem (x, y, s)
    """
    import copy
    bAx = LPb_B - (LPA_B * x2)
    UbAx = csr_matrix(spsolve(U, bAx)).T
    Qr = copy.deepcopy(Q)
    Qr.reverse()

    ordering = range(K.f)
    for i, index in enumerate(Qr):
        tmp = ordering[rank_Af - i - 1]
        ordering[rank_Af - i - 1] = ordering[index]
        ordering[index] = tmp

    # Linear dependent case
    if rank_Af < K.f:
        UbAx = bmat([[UbAx],
                    [csc_matrix((K.f - rank_Af, 1))]]).tocsc()

    x_f = UbAx[ordering]
    x = bmat([[x_f],
              [x2]])

    y = (bmat([[cfQU, y2.T]]) * LiP).T

    s = bmat([[csc_matrix((K.f, 1))],
              [s2]])

    return x, y, s
