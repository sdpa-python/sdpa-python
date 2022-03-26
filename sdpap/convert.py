#!/usr/bin/env python
"""sdpap.convert.py

Convert CLP format to SeDuMi format
This is the module of sdpap.

September 2010, Kenta KATO
"""

from .symcone import SymCone
from scipy.sparse import csc_matrix, bmat, eye
from scipy import sparse


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

    # Get size
    size_m, size_n = A.shape

    K_f = K.f
    K_l = K.l
    K_q = sum(K.q)
    K_s = sum([i ** 2 for i in K.s])

    J_f = J.f
    J_l = J.l
    J_q = sum(J.q)
    J_s = sum([i ** 2 for i in J.s])

    # ----------------------------------------
    # Make converted K.s part
    # ----------------------------------------
    if len(K.s) > 0:
        # Make index mapping from converted index to original index
        map_sdpIndex = [0] * sum([k * (k + 1) / 2 for k in K.s])
        offset = 0
        for k in K.s:
            for i in range(k):
                for j in range(i):
                    idx1 = i * k + j
                    idx2 = j * k + i
                    map_sdpIndex[i * (i + 1) / 2 + j + offset] = (idx1, idx2)
                # i == j case
                idx = i * k + i
                map_sdpIndex[i * (i + 1) / 2 + i + offset] = (idx, idx)
            offset += k * (k + 1) / 2

        # Split matrix
        c_flq = c[0:(K_f + K_l + K_q), :]
        c_s = c[(K_f + K_l + K_q):, :]
        A_flq = A[:, 0:(K_f + K_l + K_q)]
        A_s = A[:, (K_f + K_l + K_q):]

        # To make converted matrix
        convcs_row = []
        convcs_val = []
        convAs_row = []
        convAs_col = []
        convAs_val = []
        addAs_row = []
        addAs_col = []
        col_ptr = 0
        offset_row = 0
        offset_col = 0
    else:
        map_sdpIndex = None

    for sDim in K.s:
        A_block = A_s[:, col_ptr:(col_ptr + sDim ** 2)]
        c_block = c_s[col_ptr:(col_ptr + sDim ** 2)]

        # Make conv_cs
        for (row, val) in zip(c_block.nonzero()[0], c_block.data):
            i = row // sDim
            j = row % sDim
            if i >= j:
                convcs_row.append(i * (i + 1) / 2 + j + offset_col)
                convcs_val.append(val * 2 if i > j else val)

        # Make conv_As
        for (row, col, val) in zip(A_block.nonzero()[0], A_block.nonzero()[1],
                                   A_block.data):
            i = col // sDim
            j = col % sDim
            if i >= j:
                convAs_row.append(row)
                convAs_col.append(i * (i + 1) / 2 + j + offset_col)
                convAs_val.append(val * 2 if i > j else val)

        # Make add_As
        for i in range(sDim):
            for j in range(i):
                addAs_row.extend([i * sDim + j + offset_row,
                                  j * sDim + i + offset_row])
                idx = i * (i + 1) / 2 + j + offset_col
                addAs_col.extend([idx, idx])
            # i == j case
            addAs_row.append(i * sDim + i + offset_row)
            addAs_col.append(i * (i + 1) / 2 + i + offset_col)

        col_ptr += sDim ** 2
        offset_row += sDim ** 2
        offset_col += sDim * (sDim + 1) / 2

    # ----------------------------------------
    # Make converted matrix
    # ----------------------------------------
    # A2
    if K_s > 0:
        newJ_s = offset_row
        newK_s = offset_col
        sub_A = csc_matrix((convAs_val, (convAs_row, convAs_col)),
                           shape=(size_m, newK_s))
        sub_As = csc_matrix(([1.0] * len(addAs_row), (addAs_row, addAs_col)),
                            shape=(newJ_s, newK_s))
        new_A = bmat([[A_flq, sub_A]])
    else:
        newJ_s = 0
        newK_s = 0
        new_A = A

    if not sparse.isspmatrix_csc(new_A):
        new_A = new_A.tocsc()

    list_A2 = []
    list_c2 = []
    if J_f > 0:
        A_f = new_A[0:(J_f), :]
        list_A2.append(-A_f.T)
        b_f = b[0:(J_f), :]
        list_c2.append([-b_f])

    if J_l > 0:
        A_l = new_A[(J_f):(J_f + J_l), :]
        list_A2.append(-A_l.T)
        b_l = b[(J_f):(J_f + J_l), :]
        list_c2.append([-b_l])

    if K_l > 0:
        list_subAl = []
        if K_f > 0:
            list_subAl.append([csc_matrix((K_f, K_l))])
        list_subAl.append([eye(K_l, K_l, format='csc')])
        if K_q + newK_s > 0:
            list_subAl.append([csc_matrix((K_q + newK_s, K_l))])
        sub_Al = bmat(list_subAl)
        list_A2.append(-sub_Al)
        list_c2.append([csc_matrix((K_l, 1))])

    if J_q > 0:
        A_q = new_A[(J_f + J_l):(J_f + J_l + J_q), :]
        list_A2.append(-A_q.T)
        b_q = b[(J_f + J_l):(J_f + J_l + J_q), :]
        list_c2.append([-b_q])

    if K_q > 0:
        list_subAq = []
        if K_f + K_l > 0:
            list_subAq.append([csc_matrix((K_f + K_l, K_q))])
        list_subAq.append([eye(K_q, K_q, format='csc')])
        if K_s > 0:
            list_subAq.append([csc_matrix((newK_s, K_q))])
        sub_Aq = bmat(list_subAq)
        list_A2.append(-sub_Aq)
        list_c2.append([csc_matrix((K_q, 1))])

    if J_s > 0:
        A_s = new_A[(J_f + J_l + J_q):, :]
        list_A2.append(-A_s.T)
        b_s = b[(J_f + J_l + J_q):, :]
        list_c2.append([-b_s])

    if K_s > 0:
        list_subAs = []
        if K_f + K_l + K_q > 0:
            list_subAs.append([csc_matrix((K_f + K_l + K_q, K_s))])
        list_subAs.append([sub_As.T])
        sub_As = bmat(list_subAs)
        list_A2.append(-sub_As)
        list_c2.append([csc_matrix((K_s, 1))])

    A2 = bmat([list_A2])
    c2 = bmat(list_c2)

    if newK_s > 0:
        conv_cs = csc_matrix((convcs_val, (convcs_row, [0] * len(convcs_row))),
                             shape=(newK_s, 1))
        b2 = -bmat([[c_flq],
                    [conv_cs]])
    else:
        b2 = -c

    K2 = SymCone()
    K2.f = J.f
    K2.l = J.l + K.l
    K2.q = J.q + K.q
    K2.s = J.s + K.s

    J2 = SymCone()
    J2.f = A2.shape[0]

    return A2, b2, c2, K2, J2, map_sdpIndex


def clp_toEQ(A, b, c, K, J):
    """Convert from CLP format to Equality standard form (Sedumi primal format)

    THIS FUNCTION HAS NOT IMPLEMENTED YET.

    Args:
      A, b, c, K, J: CLP format

    Returns:
      Equality standard form in SeDuMi primal format, (A2, b2, c2, K2, J2)
    """
    ##################################################
    # Under construction
    ##################################################
    print('THIS FUNCTION HAS NOT IMPLEMENTED YET.')
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

def result_fromEQ(x2, y2, K, J):
    """Get result of CLP from result of converted problem by clp_toEQ().

    THIS FUNCTION HAS NOT IMPLEMENTED YET.

    Args:
      x2, y2: Result of converted problem by clp_toEQ()
      K, J  : SymCone of CLP

    Returns:
      x, y  : Result of CLP, (x, y)
    """
    ##################################################
    # Under construction
    ##################################################
    print('THIS FUNCTION HAS NOT IMPLEMENTED YET.')
    return x, y

