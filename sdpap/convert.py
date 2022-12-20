#!/usr/bin/env python
"""
Convert CLP format to SeDuMi format

These routines were originally implemented in SparseCoLO (MATLAB) [1] and are
Copyright (C) 2009 Masakazu Kojima Group
Department of Mathematical and Computing Sciences, Tokyo Institute of Technology

SparseCoLO [1] implements sparsity exploitation algorithms in [2]

Python translations of routines written for SDPAP and are
Copyright (C) 2010-2022 SDPA Project

[1] http://www.opt.c.titech.ac.jp/kojima/SparseCoLO/SparseCoLO.htm
[2] Sunyoung Kim, Masakazu Kojima, Martin Mevissen and Makoto Yamashita, "Exploiting sparsity in linear and nonlinear matrix inequalities via positive semidefinite matrix completion," Mathematical Programming, 129(1), 33–68. https://doi.org/10.1007/s10107-010-0402-6

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

September 2010: `clp_toLMI` and `result_fromLMI` written by Kenta Kato

November 2022: `clp_toEQ` and `result_fromEQ` written by Usama Muneeb. Additionally,
`clp_toLMI` code simplified and `map_sdpIndex` removed so output of `clp_toEQ`
exactly matches SparseCoLO's `CoLOtoEQform.m` (which is a wrapper on `CoLOtoLMIform.m`)
"""

from .symcone import SymCone
from scipy.sparse import csc_matrix, bmat, eye
from scipy import sparse
import numpy as np
import pdb

def clp_toLMI(A, b, c, K, J):
    """Convert from CLP format to LMI standard form (SeDuMi dual format).

    Args:
      A, b, c, K, J: CLP format

    Returns:
      A tuple (A2, b2, c2, K2, J2, map_sdpIndex). Each members are:
        LMI standard form in SeDuMi dual format, (A2, b2, c2, K2, J2)
          J2 has only an attribute 'f'
        Index mapping from converted matrix to original matrix map_sdpIndex
    """
    # Type check
    if not sparse.isspmatrix_csc(A):
        A = A.tocsc()
    if not sparse.isspmatrix_csc(b):
        b = b.tocsc()
    if not sparse.isspmatrix_csc(c):
        c = c.tocsc()
    if not K.check_validity() or not J.check_validity():
        print('clp_toLMI(): K or J is invalid.')
        return

    K_f = K.f
    K_l = K.l
    K_q = sum(K.q)
    K_s = sum([i ** 2 for i in K.s])

    J_f = J.f
    J_l = J.l
    J_q = sum(J.q)
    J_s = sum([i ** 2 for i in J.s])

    newJ_s = 0
    newK_s = 0
    new_A = A

    if not sparse.isspmatrix_csc(new_A):
        new_A = new_A.tocsc()

    list_A2 = []
    list_c2 = []
    if J_f > 0:
        print("Extracting J.f")
        A_f = new_A[0:(J_f), :]
        list_A2.append(-A_f.T)
        b_f = b[0:(J_f), :]
        list_c2.append([-b_f])

    if J_l > 0:
        print("Extracting J.l")
        A_l = new_A[(J_f):(J_f + J_l), :]
        list_A2.append(-A_l.T)
        b_l = b[(J_f):(J_f + J_l), :]
        list_c2.append([-b_l])

    if J_q > 0:
        print("Extracting J.q")
        A_q = new_A[(J_f + J_l):(J_f + J_l + J_q), :]
        list_A2.append(-A_q.T)
        b_q = b[(J_f + J_l):(J_f + J_l + J_q), :]
        list_c2.append([-b_q])

    if J_s > 0:
        print("Extracting J.s")
        A_s = new_A[(J_f + J_l + J_q):, :]
        list_A2.append(-A_s.T)
        b_s = b[(J_f + J_l + J_q):, :]
        list_c2.append([-b_s])

    # in `CoLOtoLMIform.m` (MATLAB), free cone in K is disregarded

    if K_l > 0:
        print("Extracting K.l")
        list_subAl = []
        # preceding block
        if K_f > 0:
            list_subAl.append([csc_matrix((K_f, K_l))])
        # this block
        list_subAl.append([eye(K_l, K_l, format='csc')])
        # succeeding block
        if K_q + K_s > 0:
            list_subAl.append([csc_matrix((K_q + K_s, K_l))])
        sub_Al = bmat(list_subAl) # `vstack` equivalent for sparse matrix
        list_A2.append(-sub_Al)
        list_c2.append([csc_matrix((K_l, 1))])

    if K_q > 0:
        print("Extracting K.q")
        list_subAq = []
        # preceding block
        if K_f + K_l > 0:
            list_subAq.append([csc_matrix((K_f + K_l, K_q))])
        # this block
        list_subAq.append([eye(K_q, K_q, format='csc')])
        # succeeding block
        if K_s > 0:
            list_subAq.append([csc_matrix((K_s, K_q))])
        sub_Aq = bmat(list_subAq) # `vstack` equivalent for sparse matrix
        list_A2.append(-sub_Aq)
        list_c2.append([csc_matrix((K_q, 1))])

    if K_s > 0:
        print("Extracting K.s")
        list_subAs = []
        # preceding block
        if K_f + K_l + K_q > 0:
            list_subAs.append([csc_matrix((K_f + K_l + K_q, K_s))])
        # this block
        list_subAs.append([eye(K_s, K_s, format='csc')])
        # no succeeding block required
        sub_As = bmat(list_subAs) # `vstack` equivalent for sparse matrix
        list_A2.append(-sub_As)
        list_c2.append([csc_matrix((K_s, 1))])

    A2 = bmat([list_A2], format='csr')
    c2 = bmat(list_c2) # analog of bVect = [bVectFree; bVectLP; bVectSOCP; bVectSDP];

    if K_s > 0:
        A2_flq = A2[:K_f, :]
        A2_s = A2[K_f:, :]

        size_m, size_n = A2_s.shape

        b2_flq = c[:K_f, :] # analog of `cVect`
        b2_s = c[K_f:, :] # analog of `cVectSDP`

        # `A2_s` is a block diagonal (and hence a sparse matrix)

        col_ptr = 0

        for sDim in K.s:
            # A2_s.shape = (49, 56)
            A_block = A2_s[col_ptr:(col_ptr + sDim ** 2),:]
            b_block = b2_s[col_ptr:(col_ptr + sDim ** 2)]

            LMat = np.tril(np.ones((sDim,sDim)),-1) * 2 + np.eye(sDim)
            Lvect = LMat.flatten('F')

            LMat = np.diag(Lvect)
            nzCol = ~np.all(LMat==0, axis=0)
            LMat = LMat[:, nzCol] # (49, 28)

            LMat = np.reshape(LMat, (sDim, sDim, -1))
            LMat = LMat + np.moveaxis(LMat, 0, 1)
            LMat /= 2
            LMat = np.reshape(LMat, (sDim ** 2, -1))

            addACols = np.matmul(LMat.T, A_block.todense())
            addbVect = np.matmul(LMat.T, b_block.todense())

            A2 = bmat([[A2_flq],
                    [addACols]])

            b2 = -bmat([[b2_flq],
                        [addbVect]])

    else:
        b2 = -c

    K2 = SymCone()
    # K2.f = J.f
    # K2.l = J.l + K.l
    # K2.q = J.q + K.q
    # K2.s = J.s + K.s

    J2 = SymCone()
    # J2.f = A2.shape[0]

    return A2, b2, c2, K2, J2


def clp_toEQ(A, b, c, K, J):
    """Convert from CLP format to Equality standard form (Sedumi primal format)

    Args:
      A, b, c, K, J: CLP format

    Returns:
      Equality standard form in SeDuMi primal format, (A2, b2, c2, K2, J2)
    """

    if not K.check_validity() or not J.check_validity():
        print('clp_toEQ(): K or J is invalid.')
        return

    J_f = J.f
    J_l = J.l
    J_q = sum(J.q)
    J_s = sum([i ** 2 for i in J.s])

    if (J_l + J_q + J_s) == 0:
        print("LOP to be converted into equality standard form is already equality standard form")
        return A, b, c, K, J
    else:
        A2, b2, c2, K2, J2, map_sdpIndex = clp_toLMI(-A.T, -c, -b, J, K)
        return A2, b2, c2, K2, J2


def result_fromLMI(x2, y2, K, J, map_sdpIndex):
    """Get result of CLP from result of converted problem by clp_toLMI().

    Args:
      x2, y2: Result of converted problem by clp_toLMI()
      K, J  : SymCone of CLP
      map_sdpIndex: Index mapping from converted matrix to original matrix

    Returns:
      Result of CLP, (x, y)
    """
    # Type check
    if not sparse.isspmatrix_csc(x2):
        x2 = x2.tocsc()
    if not sparse.isspmatrix_csc(y2):
        y2 = y2.tocsc()

    if len(K.s) > 0:
        K_f = K.f
        K_l = K.l
        K_q = sum(K.q)
        K_s = sum([i ** 2 for i in K.s])

        x_flq = y2[:(K.f + K.l + K_q)]
        x_s = y2[(K.f + K.l + K_q):]
        x_row = []
        x_val = []
        for (row, val) in zip(x_s.nonzero()[0], x_s.data):
            (idx1, idx2) = map_sdpIndex[row]
            if idx1 != idx2:
                x_row.extend([idx1, idx2])
                x_val.extend([val, val])
            else:
                x_row.append(idx1)
                x_val.append(val)

            x = bmat([[x_flq],
                      [csc_matrix((x_val, (x_row, [0] * len(x_row))),
                                   shape=(K_s, 1))]])
    else:
        x = y2

    size_Jq = sum(J.q)
    size_Js = sum([z ** 2 for z in J.s])
    size_Kq = sum(K.q)

    y_list = []
    if J.f + J.l > 0:
        yfl = x2[:(J.f + J.l)]
        y_list.append([yfl])

    if size_Jq > 0:
        yq = x2[(J.f + J.l + K.l):(J.f + J.l + K.l + size_Jq)]
        y_list.append([yq])

    if size_Js > 0:
        ys = x2[(J.f + J.l + K.l + size_Jq + size_Kq):
                (J.f + J.l + K.l + size_Jq + size_Kq + size_Js)]
        y_list.append([ys])

    y = bmat(y_list)
    return x, y

def result_fromEQ(x2, y2, K, J, map_sdpIndex):
    """Get result of CLP from result of converted problem by clp_toEQ().

    Args:
      x2, y2: Result of converted problem by clp_toEQ()
      K, J  : SymCone of CLP

    Returns:
      x, y  : Result of CLP, (x, y)
    """

    return x2, y2

    # return x, y

