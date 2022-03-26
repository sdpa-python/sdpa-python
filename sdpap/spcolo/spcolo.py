#!/usr/bin/env python
"""spcolo.py

A module of spcolo
Convert problem by exploiting sparsity

August 2010, Kenta KATO

December 2010: modified for scipy
"""

__all__ = ['dconv_basisrep', 'rconv_matdecomp',
           'dconv_basisresult', 'rconv_decompresult']

from . import asputils
from .clique import CliqueSet
from sdpap.symcone import SymCone
from scipy.sparse import csc_matrix, bmat, eye
from scipy import sparse

#from cvxopt import matrix, spmatrix, sparse
#import cvxopt


def dconv_basisrep(A, b, c, K, J):
    """The d-space conversion method using basis representation

    Args:
      A, b, c, K, J: CLP format

    Returns:
      A converted problem in CLP format and list of CliqueSet,
        (A2, b2, c2, K2, J2, list_cliqueSet)
    """
    if len(K.s) == 0:
        print('No conversion becsuse K.s is empty.')
        import copy
        A2 = copy.deepcopy(A)
        b2 = copy.deepcopy(b)
        c2 = copy.deepcopy(c)
        K2 = K.SymCone()
        J2 = J.SymCone()
        return A2, b2, c2, K2, J2, None

    sDimTotal = sum([x ** 2 for x in K.s])
    size_m, size_n = A.shape
    numSDPcones = len(K.s)

    list_cliqueSet = []

    # When a block is not converted
    list_addcs = []
    list_addAs = []
    list_addKs = []

    # Additional J.s part
    list_addAclq = []
    list_addJs = []

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
    # Get maximal cliques and make additional submatrix
    # --------------------------------------------------
    col_ptr = K.f + K.l + sum(K.q)
    offset_row = 0
    offset_col = 0
    # Additional c_f
    addcf_row = []
    addcf_val = []
    # Additional A_f
    addAf_row = []
    addAf_col = []
    addAf_val = []
    # Additional J.s part of A_f
    addAclq_row = []
    addAclq_col = []

    for sDim in K.s:
        A_block = A[:, col_ptr:(col_ptr + sDim ** 2)]
        c_block = c[col_ptr:(col_ptr + sDim ** 2)]

        # Get Aggregate Sparsity Pattern
        asp = asputils.getASP(A_block, c_block, (sDim,))
        # print asp

        # Get maximal clique and add to clique list
        cliqueSet = asputils.cliques_fromASP(asp)
        if cliqueSet == None:
            # Doesn't need to convert this block
            list_addcs.append([c_block])
            list_addAs.append(A_block)
            list_addKs.append(sDim)
            list_cliqueSet.append(None)
            continue

        list_cliqueSet.append(cliqueSet)
        # print cliqueSet.map_convIndex

        # ----------------------------------------
        # Make additional submatrix
        # ----------------------------------------

        # Make addAclq
        add_Js = [len(clq) for clq in cliqueSet.cliques]
        addJs_size = sum([x ** 2 for x in add_Js])
        addAclq_row.extend([i + offset_row for i in range(addJs_size)])
        for clq in cliqueSet.cliques:
            addAclq_col.extend([cliqueSet.map_convIndex[i*sDim+j] + offset_col
                                for i in clq for j in clq])

        # Make addcf
        for (row, val) in zip(c_block.nonzero()[0], c_block.data):
            i = row // sDim
            j = row % sDim
            if i >= j:
                addcf_row.append(cliqueSet.map_convIndex[row] + offset_col)
                addcf_val.append(val * 2 if i > j else val)

        #print addcf_row

        # Make addAf
        for (row, col, val) in zip(A_block.nonzero()[0], A_block.nonzero()[1],
                                   A_block.data):
            i = col // sDim
            j = col % sDim
            if i >= j:
                addAf_row.append(row)
                addAf_col.append(cliqueSet.map_convIndex[col] + offset_col)
                addAf_val.append(val * 2 if i > j else val)


        list_addJs.extend(add_Js)

        col_ptr += sDim ** 2
        offset_row += addJs_size
        offset_col += cliqueSet.num_element

    # --------------------------------------------------
    # Make converted matrix
    # --------------------------------------------------
    size_addJs = offset_row
    size_addKf = offset_col
    size_addKs = sum(list_addKs)

    list_A = []
    list_c = []

    # J.f, J.l, J.q, J.s part
    # K.f part
    if K.f > 0:
        A_f = A[:, 0:K.f]
        list_A.append(A_f)
        c_f = c[0:K.f]
        list_c.append([c_f])

    # Additional K.f part
    if len(addAf_val) > 0:
        add_Af = csc_matrix((addAf_val, (addAf_row, addAf_col)),
                            shape=(size_m, size_addKf))
        add_cf = csc_matrix((addcf_val, (addcf_row, [0] * len(addcf_row))),
                            shape=(size_addKf, 1))
        list_A.append(add_Af)
        list_c.append([add_cf])

    # K.l, K.q part
    if K.l > 0 or len(K.q) > 0:
        A_lq = A[:, K.f:(K.f + K.l + sum(K.q))]
        list_A.append(A_lq)
        c_lq = c[K.f:(K.f + K.l + sum(K.q))]
        list_c.append([c_lq])

    # Additional K.s part
    if len(addAf_val) > 0:
        list_A.extend(list_addAs)
        list_c.extend(list_addcs)

    # Concatenate matrix
    A2 = bmat([list_A])
    c2 = bmat(list_c)

    # Additional J.s part
    if len(addAf_val) > 0:
        list_A = []
        add_Aclq = csc_matrix(([1.0] * len(addAclq_row),
                               (addAclq_row, addAclq_col)),
                              shape=(size_addJs, size_addKf))
        if K.f > 0:
            list_A.append(csc_matrix((size_addJs, K.f)))

        list_A.append(add_Aclq)

        if K.l > 0 or len(K.q) > 0:
            list_A.append(csc_matrix((size_addJs, K.l + sum(K.q))))

        if len(list_addAs) > 0:
            list_A.append(csc_matrix((size_addJs, size_addKs)))

        add_A = bmat([list_A])
        A2 = bmat([[A2], [add_A]])
        b2 = bmat([[b], [csc_matrix((size_addJs, 1))]])
    else:
        b2 = b

    # --------------------------------------------------
    # Make converted SymCone
    # --------------------------------------------------
    K2 = SymCone(f=K.f+size_addKf, l=K.l, q=K.q, s=tuple(list_addKs))
    J2 = SymCone(f=J.f, l=J.l, q=J.q, s=tuple(list(J.s)+list_addJs))

    return A2, b2, c2, K2, J2, list_cliqueSet


def rconv_matdecomp(A, b, c, K, J):
    """The r-space conversion method using matrix decomposition.

    Args:
      A, b, c, K, J: CoLO format

    Returns:
      A converted problem in CoLO format and list of CliqueSet,
        (A2, b2, c2, K2, J2, list_cliqueSet)
    """
    dom_A, dom_b, dom_c, dom_K, dom_J, list_cliqueSet = \
           dconv_basisrep(-A.T, -c, -b, J, K)
    A2 = -dom_A.T
    c2 = -dom_b
    b2 = -dom_c
    J2 = dom_K
    K2 = dom_J

    return A2, b2, c2, K2, J2, list_cliqueSet


def dconv_basisresult(x, y, K, J, dom_K, cliqueD):
    """Retrieve an optimal solution which solved by dconv_basisrep()

    Args:
      x, y: Results of the converted problem
      K, J: A SymCone of the original problem
      dom_K: A SymCone of the converted problem
      cliqueD: A list of CliqueSet

    Returns:
      x0, y0: Results of the original problem
    """
    # --------------------------------------------------
    # Type check
    # --------------------------------------------------
    if not sparse.isspmatrix_csc(x):
        x = x.tocsc()
    if not sparse.isspmatrix_csc(y):
        y = y.tocsc()

    x0_row = []
    x0_val = []
    if K.f > 0:
        x_f = x[:K.f]
        x0_row.extend(list(x_f.nonzero()[0]))
        x0_val.extend(list(x_f.data))

    add_xf = x[K.f:dom_K.f]

    if K.l > 0 or len(K.q) > 0:
        x_lq = x[dom_K.f:(dom_K.f + dom_K.l + sum(dom_K.q))]
        x0_row.extend(list(x_lq.nonzero()[0]))
        x0_val.extend(list(x_lq.data))

    if len(dom_K.s) > 0:
        add_xs = x[(dom_K.f + dom_K.l + sum(dom_K.q)):]

    offset_x0 = K.f + K.l + sum(K.q)
    offset_addxf = 0
    offset_addxs = 0
    index_addKs = 0
    for cliqueSet in cliqueD:
        if cliqueSet != None:
            size_clique = cliqueSet.num_element
            addxf_block = add_xf[offset_addxf:(offset_addxf + size_clique)]
            for (row, val) in zip(addxf_block.nonzero()[0], addxf_block.data):
                orig_idx = cliqueSet.map_origIndex[row]
                if orig_idx[0] != orig_idx[1]:
                    x0_row.extend([offset_x0 + orig_idx[0],
                                   offset_x0 + orig_idx[1]])
                    x0_val.extend([val, val])
                else:
                    x0_row.append(offset_x0 + orig_idx[0])
                    x0_val.append(val)

            offset_addxf += size_clique
            offset_x0 += size_clique
        else:
            size_Ks = dom_K.s[index_addKs] ** 2
            addxs_block = add_xs[offset_addxs:(offset_addxs + size_Ks)]
            addxs_row = list(addxs_block.I)
            for row in addxs_row:
                x0_row.append(offset_x0 + row)

            x0_val.extend(list(addxs_block.V))
            offset_addxs += size_Ks
            offset_x0 += size_Ks

    totalsize_x = K.f + K.l + sum(K.q) + sum([i ** 2 for i in K.s])
    totalsize_y = J.f + J.l + sum(J.q) + sum([i ** 2 for i in J.s])
    x0 = csc_matrix((x0_val, (x0_row, [0] * len(x0_row))),
                           shape=(totalsize_x, 1))
    y0 = y[:totalsize_y]

    return x0, y0


def rconv_decompresult(x, y, K, J, ran_J, cliqueR):
    """Retrieve an optimal solution which solved by rconv_matdecomp()

    Args:
      x, y: Results of the converted problem
      K, J: A SymCone of the original problem
      ran_J: A SymCone of the converted problem
      cliqueR: A list of CliqueSet

    Returns:
      x0, y0: Results of the original problem
    """
    y0, x0 = dconv_basisresult(y, x, J, K, ran_J, cliqueR)
    return x0, y0
