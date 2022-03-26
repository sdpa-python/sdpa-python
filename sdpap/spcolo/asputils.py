#!/usr/bin/env python
"""asputils.py

A module of spcolo
Subprograms using Aggregate Sparsity Pattern(ASP)

September 2010, Kenta KATO

December, 2010: modified for scipy
"""

__all__ = ['getASP', 'cliques_fromASP']

from sdpap.matdata import MatData
from . import spcoloext
from .clique import CliqueSet
from scipy.sparse import csc_matrix, tril
from scipy import array, sparse

def getASP(A, c, K_s):
    """Get Aggregate Sparsity Pattern

    Args;
      A, c: submatrix or subvector of CoLO
      K_s: A tuple of block struct of A or c

    Returns:
      Aggregate Sparsity Pattern Matrix
    """
    if not isinstance(K_s, tuple):
        raise TypeError('K_s must be tuple.')

    if len(K_s) == 0:
        print("getASP(): K_s is not assigned.")
        return

    size_asp = sum(K_s)
    sDim = sum([x ** 2 for x in K_s])
    size_m, size_n = A.shape

    asp_I = set(c.nonzero()[0])
    asp_I.update(set(A.nonzero()[1]))

    list_I = sorted(list(asp_I))
    mat_offset = 0
    vec_offset = 0
    asp_row = []
    asp_col = []
    for k in K_s:
        block_I = list_I[vec_offset:(vec_offset + k ** 2)]
        asp_row.extend([((i - vec_offset) // k) + mat_offset
                        for i in block_I])
        asp_col.extend([((i - vec_offset) % k) + mat_offset
                        for i in block_I])
        mat_offset += k
        vec_offset += k ** 2

    ret_mat = csc_matrix(([1] * len(asp_row), (asp_row, asp_col)),
                         shape=(size_asp, size_asp))

    return ret_mat


def symb_cholesky(asp):
    """Symbolic Cholesky Factorization

    Args:
      asp: Aggregated Sparsity Pattern Matrix

    Returns:
      Indices and indptr of Lower triangular matrix L (indL, ptrL)
    """
    A = tril(asp, k=-1, format='csc')
    A.sort_indices()
    parent = dict()
    indA = A.indices
    ptrA = A.indptr
    setL = []
    for i in range(A.shape[0]):
        setL.append(set(indA[ptrA[i]:ptrA[i+1]]))
        for j in parent.keys():
            if parent[j] == i:
                setL[i].update(setL[j])

        setL[i] -= set(range(i+1))
        if setL[i]:
            parent[i] = min(setL[i])

    indL = []
    ptrL = [0] * (A.shape[0] + 1)
    for i in range(A.shape[0]):
        indL.append(i)
        indL.extend(sorted(list(setL[i])))
        ptrL[i+1] = len(indL)

    #L = csc_matrix((array([1.0] * len(indL)),
    #                array(indL),
    #                array(ptrL)),
    #               shape=A.shape)
    return indL, ptrL


def cliques_fromASP(asp):
    """Get cliques from Aggregate Sparsity Pattern Matrix

    Args:
      asp: Aggregate Sparsity Pattern Matrix

    Returns:
      orig_cliques: A list of list. Each list has indexes of clique.
    """
    from scipy import ix_

    size_n = asp.shape[0]
    A = asp + sparse.eye(size_n, size_n, format='csc')
    #A = asp + (2 * size_n + 1) * sparse.eye(size_n, size_n, format='csc')
    data_A = MatData(mat=A.tocsc().sorted_indices())

    # Minimum degree ordering
    # using spooles
    ordering = spcoloext.ordering_mmd(data_A)
    if len(ordering) == 0:
        return None

    if A.shape != (1, 1):
        #data_A = MatData(mat=A.tocsc()[ix_(ordering, ordering)])
        A = A.tocsc()[ix_(ordering, ordering)]
    else:
        #data_A = MatData(mat=A.tocsc())
        A = A.tocsc()

    # Cholesky factorization
    #L = symb_cholesky(A)
    indL, ptrL = symb_cholesky(A)
    #values, rowind, colptr = spcoloext.cholesky(data_A)
    #data_R = MatData(values=values, rowind=rowind, colptr=colptr,
    #                 size=(size_n,size_n))

    #R = data_R.tomatrix()
    #L = R.T

    # Finding the maximal cliques
    #clique_set = [set(L[:,col].nonzero()[0]) for col in range(size_n)]
    clique_set = [set(indL[ptrL[col]:ptrL[col+1]]) for col in range(size_n)]
    maxclique_idx = [0]
    for i in range(1, size_n):
        for j in range(i):
            if clique_set[i] <= clique_set[j]:
                break
        else:
            maxclique_idx.append(i)

    cliqueSet = CliqueSet(size_n)
    for i in maxclique_idx:
        clq = sorted([ordering[j] for j in clique_set[i]])
        cliqueSet.append_clique(clq)

    # print cliqueSet.cliques

    return cliqueSet
