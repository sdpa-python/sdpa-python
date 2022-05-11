#!/usr/bin/env python
"""fileio.py

Input or output SDPA sparse format (dat-s) file
This is the module of sdpap.

September 2010, Kenta KATO

December 2010: modified for scipy
"""

__all__ = ['readproblem', 'writeproblem', 'fromsdpa', 'tosdpa']

from . import convert
from .symcone import SymCone
from scipy.sparse import csc_matrix, csr_matrix
from scipy import sparse
from numpy import matrix


def readproblem(filename):
    """Read problem file and make CLP

    Args:
      filename: A string of file pass for CLP format problem
        Input file format:
          J.f
          J.l
          J.q        -> Each numbers are split by " "(white space).
                          If J.q = (), write 0.
          J.s        -> Each numbers are split by " ". if J.s = (), write 0.
          K.f
          K.l
          K.q        -> Each numbers are split by " ". if K.q = (), write 0.
          K.s        -> Each numbers are split by " ". if K.s = (), write 0.
          ROW COL row col val  -> Each numbers are split by " ".
          ROW COL row col val
           :   :   :   :   :
           :   :   :   :   :

          (ROW,COL,row,col,val) means:
          The value of (row,col) element in (ROW,COL) block of A is val
            (ROW,COL,row,col are indexed 1,2,3,...)
          If ROW = 0: the value of col's element in
                        COL's block of c is val (row = 0(recommended))
          If COL = 0: the value of row's element in
                        ROW's block of b is val (col = 0(recommended))

    Returns:
      A tuple of CLP problem input (A, b, c, K, J)
    """
    if not isinstance(filename, str):
        print("readproblem(): filename must be string.")
        return

    fp = open(filename, "r")
    # Read J
    J_f = int(fp.readline())
    J_l = int(fp.readline())
    line = [int(x) for x in fp.readline().rstrip().split(" ")]
    J_q = tuple(line) if line[0] != 0 else ()
    line = [int(x) for x in fp.readline().rstrip().split(" ")]
    J_s = tuple(line) if line[0] != 0 else ()
    J = SymCone(J_f, J_l, J_q, J_s)

    # Read K
    K_f = int(fp.readline())
    K_l = int(fp.readline())
    line = [int(x) for x in fp.readline().rstrip().split(" ")]
    K_q = tuple(line) if line[0] != 0 else ()
    line = [int(x) for x in fp.readline().rstrip().split(" ")]
    K_s = tuple(line) if line[0] != 0 else ()
    K = SymCone(K_f, K_l, K_q, K_s)

    start_row = [0]
    if J.f != 0:
        start_row.append(J.f)

    if J.l != 0:
        start_row.append(start_row[-1] + J.l)

    for k in J.q:
        start_row.append(start_row[-1] + k)

    for k in J.s:
        start_row.append(start_row[-1] + k ** 2)

    start_col = [0]
    if K.f != 0:
        start_col.append(K.f)

    if K.l != 0:
        start_col.append(start_col[-1] + K.l)

    for k in K.q:
        start_col.append(start_col[-1] + k)

    for k in K.s:
        start_col.append(start_col[-1] + k ** 2)

    size_row = J.f + J.l + sum(J.q) + sum([k ** 2 for k in J.s])
    size_col = K.f + K.l + sum(K.q) + sum([k ** 2 for k in K.s])

    # Read I,J,i,j,v
    A_row = []
    A_col = []
    A_val = []
    b_index = []
    b_val = []
    c_index = []
    c_val = []

    for line in fp.readlines():
        ROW, COL, row, col, val = line.rstrip().split(" ")[0:5]
        if int(ROW) == 0:
            c_index.append(start_col[int(COL) - 1] + int(col) - 1)
            c_val.append(float(val))
        elif int(COL) == 0:
            b_index.append(start_row[int(ROW) - 1] + int(row) - 1)
            b_val.append(float(val))
        else:
            A_row.append(start_row[int(ROW) - 1] + int(row) - 1)
            A_col.append(start_col[int(COL) - 1] + int(col) - 1)
            A_val.append(float(val))

    A = csc_matrix((A_val, (A_row, A_col)), shape=(size_row, size_col))
    b = csc_matrix((b_val, (b_index, [0] * len(b_index))), shape=(size_row, 1))
    c = csc_matrix((c_val, (c_index, [0] * len(c_index))), shape=(size_col, 1))

    fp.close()

    return A, b, c, K, J


def writeproblem(filename, A, b, c, K, J, accuracy="%+8.16e"):
    """Write CLP to file

    Args:
      filename: A string of file pass
      A, b, c: Scipy matrices to denote the CLP
      K, J: Symcone object to denote the CLP
      accuracy: Print format
    """
    if not isinstance(filename, str):
        raise ValueError("writeproblem(): filename must be string.")

    if not sparse.isspmatrix_csc(b):
        b = csc_matrix(b)

    if not sparse.isspmatrix_csc(c):
        c = csc_matrix(c)

    if not sparse.isspmatrix_csc(A):
        A = csc_matrix(A)

    fp = open(filename, "w")
    # Write J
    fp.write(str(J.f) + "\n")
    fp.write(str(J.l) + "\n")
    if len(J.q) > 0:
        fp.write(" ".join([str(num) for num in J.q]) + "\n")
    else:
        fp.write("0\n")
    if len(J.s) > 0:
        fp.write(" ".join([str(num) for num in J.s]) + "\n")
    else:
        fp.write("0\n")

    # Write K
    fp.write(str(K.f) + "\n")
    fp.write(str(K.l) + "\n")
    if len(K.q) > 0:
        fp.write(" ".join([str(num) for num in K.q]) + "\n")
    else:
        fp.write("0\n")
    if len(K.s) > 0:
        fp.write(" ".join([str(num) for num in K.s]) + "\n")
    else:
        fp.write("0\n")

    start_row = [0]
    if J.f != 0:
        start_row.append(J.f)

    if J.l != 0:
        start_row.append(start_row[-1] + J.l)

    for k in J.q:
        start_row.append(start_row[-1] + k)

    for k in J.s:
        start_row.append(start_row[-1] + k ** 2)

    start_col = [0]
    if K.f != 0:
        start_col.append(K.f)

    if K.l != 0:
        start_col.append(start_col[-1] + K.l)

    for k in K.q:
        start_col.append(start_col[-1] + k)

    for k in K.s:
        start_col.append(start_col[-1] + k ** 2)

    size_row = J.f + J.l + sum(J.q) + sum([k ** 2 for k in J.s])
    size_col = K.f + K.l + sum(K.q) + sum([k ** 2 for k in K.s])

    # Split matrix
    c_block = [c[start_col[i]:start_col[i+1]]
               for i in range(len(start_col) - 1)]
    b_block = [b[start_row[i]:start_row[i+1]]
               for i in range(len(start_row) - 1)]
    tmp_block = [A[start_row[i]:start_row[i+1], :]
                 for i in range(len(start_row) - 1)]

    A_block = []
    for T in tmp_block:
        A_block.append([T[:, start_col[i]:start_col[i+1]]
                        for i in range(len(start_col) - 1)])

    # Write c
    for COL in range(len(start_col) - 1):
        block = c_block[COL]
        list_row = ["0"] * len(list(block.nonzero()[0]))
        list_col = [str(x + 1) for x in list(block.nonzero()[0])]
        list_val = [accuracy % x for x in list(block.data)]
        for (row, col, val) in zip(list_row, list_col, list_val):
            fp.write(" ".join(["0", str(COL + 1), row, col, val]) + "\n")

    # Write b
    for ROW in range(len(start_row) - 1):
        block = b_block[ROW]
        list_row = [str(x + 1) for x in list(block.nonzero()[0])]
        list_col = ["0"] * len(list(block.nonzero()[0]))
        list_val = [accuracy % x for x in list(block.data)]
        for (row, col, val) in zip(list_row, list_col, list_val):
            fp.write(" ".join([str(ROW + 1), "0", row, col, val]) + "\n")

    # Write A
    for ROW in range(len(start_row) - 1):
        for COL in range(len(start_col) - 1):
            block = A_block[ROW][COL]
            list_row = [str(x + 1) for x in list(block.nonzero()[0])]
            list_col = [str(x + 1) for x in list(block.nonzero()[1])]
            list_val = [accuracy % x for x in list(block.data)]
            for (row, col, val) in zip(list_row, list_col, list_val):
                fp.write(" ".join([str(ROW + 1), str(COL + 1),
                                   row, col, val]) + "\n")

    fp.close()
    return


def fromsdpa(filename):
    """Convert from SDPA sparse format to CLP format

    Args:
      filename: A string of file pass for SDPA sparse format

    Returns:
      A tuple of CLP problem input, (A, b, c, K, J)
    """
    if not isinstance(filename, str):
        print('fromsdpa(): filename must be string')
        return

    K = SymCone()
    J = SymCone()

    # skip comment
    fp = open(filename, 'r')
    line = fp.readline()
    while line[0] == '*' or line[0] == '"':
        line = fp.readline()

    # read mDim
    J.f = int(line.strip().split(' ')[0])

    # read nBlock
    line = fp.readline()
    nBlock = int(line.strip().split(' ')[0])

    # read blockStruct
    blockStruct = []
    structSOCP = []
    structSDP = []
    line = fp.readline()
    blockStruct = line.strip().split(' ')
    while blockStruct.count('') > 0:
        blockStruct.remove('')

    blockSize = []
    blockType = []
    for i in range(nBlock):
        item = blockStruct[i]
        size = abs(int(item[0:-1])) if item[-1] in 'FLQS' else abs(int(item))
        blockSize.append(size)

        if item[-1] == 'F':
            K.f += size
            blockType.append('F')
        elif item[-1] == 'L':
            K.l += size
            blockType.append('L')
        elif item[-1] == 'Q':
            structSOCP.append(size)
            blockType.append('Q')
        elif item[-1] == 'S':
            structSDP.append(size)
            blockType.append('S')
        else:
            st = int(item)
            if st > 0:
                structSDP.append(st)
                blockType.append('S')
            else:
                K.l -= st
                blockType.append('L')

    K.q = tuple(structSOCP)
    K.s = tuple(structSDP)

    # read b
    line = fp.readline()
    line = line.strip()
    line = line.strip('{}()')
    if ',' in line:
        b_str = line.strip().split(',')
    else:
        b_str = line.strip().split()
    while b_str.count('') > 0:
        b_str.remove('')

    b = csr_matrix(matrix(list(map(float, [s for s in b_str])))).T

    # read c and A
    blockElements_c = [[] for i in range(nBlock)]
    blockElements_A = [[] for i in range(nBlock)]
    lineList = fp.readlines()
    for line in lineList:
        row, block, colI, colJ, val = line.split()[0:5]
        row = int(row.strip(',')) - 1
        block = int(block.strip(',')) - 1
        colI = int(colI.strip(',')) - 1
        colJ = int(colJ.strip(',')) - 1
        val = float(val.strip(','))
        col = colI * blockSize[block] + colJ \
              if blockType[block] == 'S' else colI
        col2 = colJ * blockSize[block] + colI \
               if blockType[block] == 'S' else colJ

        if row == -1:
            blockElements_c[block].append((col, val))
            if blockType[block] == 'S' and colI != colJ:
                blockElements_c[block].append((col2, val))
        else:
            blockElements_A[block].append((row, col, val))
            if blockType[block] == 'S' and colI != colJ:
                blockElements_A[block].append((row, col2, val))

    for block in range(nBlock):
        blockElements_c[block].sort()
        blockElements_A[block].sort(key = lambda x: x[1])

    elementsK_c = {'F':[[], []], 'L':[[], []], 'Q':[[], []], 'S':[[], []]}
    elementsK_A = {'F':[[], [], []], 'L':[[], [], []],
                   'Q':[[], [], []], 'S':[[], [], []]}
    offsetK = {'F':0, 'L':0, 'Q':0, 'S':0}

    for block in range(nBlock):
        for item in blockElements_c[block]:
            elementsK_c[blockType[block]][0].append(item[0] +
                                                    offsetK[blockType[block]])
            elementsK_c[blockType[block]][1].append(item[1])

        for item in blockElements_A[block]:
            elementsK_A[blockType[block]][0].append(item[0])
            elementsK_A[blockType[block]][1].append(item[1] +
                                                    offsetK[blockType[block]])
            elementsK_A[blockType[block]][2].append(item[2])

        offsetK[blockType[block]] += blockSize[block] ** 2 \
                                     if blockType[block] == 'S' \
                                     else blockSize[block]

    elements_c = [[], []]
    elements_A = [[], [], []]
    offset = 0

    if K.f > 0:
        elements_c[0].extend(elementsK_c['F'][0])
        elements_c[1].extend(elementsK_c['F'][1])
        elements_A[0].extend(elementsK_A['F'][0])
        elements_A[1].extend(elementsK_A['F'][1])
        elements_A[2].extend(elementsK_A['F'][2])
        offset += K.f

    if K.l > 0:
        elements_c[0].extend([x + offset for x in elementsK_c['L'][0]])
        elements_c[1].extend(elementsK_c['L'][1])
        elements_A[0].extend(elementsK_A['L'][0])
        elements_A[1].extend([x + offset for x in elementsK_A['L'][1]])
        elements_A[2].extend(elementsK_A['L'][2])
        offset += K.l

    if True: #K.q > 0:
        elements_c[0].extend([x + offset for x in elementsK_c['Q'][0]])
        elements_c[1].extend(elementsK_c['Q'][1])
        elements_A[0].extend(elementsK_A['Q'][0])
        elements_A[1].extend([x + offset for x in elementsK_A['Q'][1]])
        elements_A[2].extend(elementsK_A['Q'][2])
        offset += sum(K.q)

    if True: #K.s > 0:
        elements_c[0].extend([x + offset for x in elementsK_c['S'][0]])
        elements_c[1].extend(elementsK_c['S'][1])
        elements_A[0].extend(elementsK_A['S'][0])
        elements_A[1].extend([x + offset for x in elementsK_A['S'][1]])
        elements_A[2].extend(elementsK_A['S'][2])
        offset += sum([x ** 2 for x in K.s])

    size_K = offset

    c = csc_matrix((elements_c[1], (elements_c[0], [0] * len(elements_c[0]))),
                   shape=(size_K, 1))
    A = csc_matrix((elements_A[2], (elements_A[0], elements_A[1])),
                   shape=(J.f, size_K))

    return A, b, c, K, J


def tosdpa(filename, A, b, c, K, J, accuracy="%+8.16e"):
    """Convert from SeDuMi format to SDPA sparse format

    If J.l or J.q or J.s are not empty, CLP format will be converted
      to SeDuMi format first.

    Args:
      filename: A string of file pass
      A, b, c: Scipy matrices to denote the CLP
      K, J: Symcone object to denote the CLP
      accuracy: Print format (e.g. '%8.16e')
    """
    if accuracy[0] != "%":
        raise ValueError("accuracy must start with %% (e.g. %8.16e) \n")

    if not isinstance(K, SymCone) or not isinstance(J, SymCone):
        raise ValueError('K and J must be an instance of SymCone.')

    if not K.check_validity() or not J.check_validity():
        return

    if J.l > 0 or len(J.q) > 0 or len(J.s) > 0:
        raise ValueError("SymCone J must only have attribute 'f'")

    if not sparse.isspmatrix_csc(A):
        A = A.tocsc()
    if not sparse.isspmatrix_csc(b):
        b = b.tocsc()
    if not sparse.isspmatrix_csc(c):
        c = c.tocsc()

    # Note that primal-dual is reverse in SeDuMi.
    # So c2 must be -c2.
    # In addition, A2 should be passed in the transposed style.
    c2 = -c
    b2 = -b
    A2 = -A

    fp = open(filename, "w")

    size_row, size_col = A2.shape
    # write mDim
    fp.write(str(size_row) + "\n")

    # write nBlock
    fp.write(str(int(K.f > 0) + int(K.l > 0) + len(K.q) + len(K.s)) + "\n")

    print(J)
    print(K)

    # write blockStruct
    blockStruct = []
    if K.f > 0:
        blockStruct.append(str(K.f) + "F")

    if K.l > 0:
        blockStruct.append(str(K.l) + "L")

    if len(K.q) > 0:
        blockStruct.extend([str(x) + "Q" for x in K.q])

    if len(K.s) > 0:
        blockStruct.extend([str(x) + "S" for x in K.s])

    fp.write(" ".join(blockStruct) + "\n")

    # write b
    fp.write(" ".join([accuracy % b2[i,0] for i in range(b2.shape[0])]) + "\n")

    len_q = sum(K.q)
    len_s = sum([x ** 2 for x in K.s])
    # write c
    if K.f > 0:
        c_f = c2[0:K.f, :]
    if K.l > 0:
        c_l = c2[K.f:(K.f + K.l), :]
    if len_q > 0:
        c_q = c2[(K.f + K.l):(K.f + K.l + len_q), :]
    if len_s > 0:
        c_s = c2[(K.f + K.l + len_q):, :]

    setBlock = 1
    if K.f > 0:
        list_row = [str(x + 1) for x in list(c_f.nonzero()[0])]
        list_val = [accuracy % x for x in list(c_f.data)]
        length = len(list_row)
        for i in range(length):
            fp.write(" ".join(("0", str(setBlock), list_row[i],
                               list_row[i], list_val[i])) + "\n")

        setBlock += 1

    if K.l > 0:
        list_row = [str(x + 1) for x in list(c_l.nonzero()[0])]
        list_val = [accuracy % x for x in list(c_l.data)]
        length = len(list_row)
        for i in range(length):
            fp.write(" ".join(("0", str(setBlock), list_row[i],
                               list_row[i], list_val[i])) + "\n")

        setBlock += 1

    if len_q > 0:
        offset = 0
        for blockSize in K.q:
            list_row = [str(x + 1) for x in
                        list(c_q[offset:(offset + blockSize), :].nonzero()[0])]
            list_val = [accuracy % x for x in
                        list(c_q[offset:(offset + blockSize), :].data)]
            length = len(list_row)
            for i in range(length):
                fp.write(" ".join(("0", str(setBlock), list_row[i],
                                   list_row[i], list_val[i])) + "\n")

            setBlock += 1
            offset += blockSize

    if len_s > 0:
        offset = 0
        for blockSize in K.s:
            list_row = list(c_s[offset:(offset + blockSize * blockSize),
                                :].nonzero()[0])
            list_val = [accuracy % x for x in
                        list(c_s[offset:(offset + blockSize * blockSize),
                                 :].data)]
            length = len(list_row)
            for i in range(length):
                setCol_row = (list_row[i] // blockSize) + 1
                setCol_col = (list_row[i] % blockSize) + 1
                if setCol_row <= setCol_col:
                    fp.write(" ".join(("0", str(setBlock), str(setCol_row),
                                       str(setCol_col), list_val[i])) + "\n")

            setBlock += 1
            offset += blockSize * blockSize

    # write A
    if K.f > 0:
        A_f = A2[:, 0:K.f]
    if K.l > 0:
        A_l = A2[:, K.f:(K.f + K.l)]
    if len_q > 0:
        A_q = A2[:, (K.f + K.l):(K.f + K.l + len_q)]
    if len_s > 0:
        A_s = A2[:, (K.f + K.l + len_q):]

    setBlock = 1
    if K.f > 0:
        list_row = [str(x + 1) for x in list(A_f.nonzero()[0])]
        list_col = [str(x + 1) for x in list(A_f.nonzero()[1])]
        list_val = [accuracy % x for x in list(A_f.data)]
        length = len(list_row)
        for i in range(length):
            fp.write(" ".join((list_row[i], str(setBlock), list_col[i],
                               list_col[i], list_val[i])) + "\n")

        setBlock += 1

    if K.l > 0:
        list_row = [str(x + 1) for x in list(A_l.nonzero()[0])]
        list_col = [str(x + 1) for x in list(A_l.nonzero()[1])]
        list_val = [accuracy % x for x in list(A_l.data)]
        length = len(list_row)
        for i in range(length):
            fp.write(" ".join((list_row[i], str(setBlock), list_col[i],
                               list_col[i], list_val[i])) + "\n")

        setBlock += 1

    if len_q > 0:
        offset = 0
        for blockSize in K.q:
            list_row = [str(x + 1) for x in
                        list(A_q[:, offset:(offset + blockSize)].nonzero()[0])]
            list_col = [str(x + 1) for x in
                        list(A_q[:, offset:(offset + blockSize)].nonzero()[1])]
            list_val = [accuracy % x for x in
                        list(A_q[:, offset:(offset + blockSize)].data)]
            length = len(list_row)
            for i in range(length):
                fp.write(" ".join((list_row[i], str(setBlock), list_col[i],
                                   list_col[i], list_val[i])) + "\n")

            setBlock += 1
            offset += blockSize

    if len_s > 0:
        offset = 0
        for blockSize in K.s:
            list_row = [str(x + 1) for x in
                        list(A_s[:, offset:(offset + blockSize * blockSize)].
                             nonzero()[0])]
            list_col = list(A_s[:, offset:(offset + blockSize * blockSize)].
                            nonzero()[1])
            list_val = [accuracy % x for x in
                        list(A_s[:, offset:(offset + blockSize * blockSize)].
                             data)]
            length = len(list_row)
            for i in range(length):
                setCol_row = (list_col[i] // blockSize) + 1
                setCol_col = (list_col[i] % blockSize) + 1
                if setCol_row <= setCol_col:
                    fp.write(" ".join((list_row[i], str(setBlock),
                                       str(setCol_row), str(setCol_col),
                                       list_val[i])) + "\n")

            setBlock += 1
            offset += blockSize * blockSize

    fp.close()
    return


# ==================================================
# Functions to write result file
# ==================================================
def write_version(fp):
    """Write SDPAP title and version to file

    Args:
      fp: File pointer
    """
    fp.write("==================================================\n" +
             " SDPAP: SDPA Python Interface\n" +
             " SDPA group 2010-2011\n" +
             "==================================================\n")
    return

def write_parameter(fp, option):
    """Write parameters to file

    Args:
      fp: File pointer
      option: Parameters
    """
    fp.write("----- Parameters -----\n" +
             'maxIteration: ' + str(option['maxIteration']) + "\n" +
             ' epsilonStar: ' + str(option['epsilonStar']) + "\n" +
             '  lambdaStar: ' + str(option['lambdaStar']) + "\n" +
             '   omegaStar: ' + str(option['omegaStar']) + "\n" +
             '  lowerBound: ' + str(option['lowerBound']) + "\n" +
             '  upperBound: ' + str(option['upperBound']) + "\n" +
             '    betaStar: ' + str(option['betaStar']) + "\n" +
             '     betaBar: ' + str(option['betaBar']) + "\n" +
             '   gammaStar: ' + str(option['gammaStar']) + "\n" +
             ' epsilonDash: ' + str(option['epsilonDash']) + "\n" +
             ' isSymmetric: ' + str(option['isSymmetric']) + "\n" +
             '    isDimacs: ' + str(option['isDimacs']) + "\n" +
             '      xPrint: ' + str(option['xPrint']) + "\n" +
             '      yPrint: ' + str(option['yPrint']) + "\n" +
             '      sPrint: ' + str(option['sPrint']) + "\n" +
             '    infPrint: ' + str(option['infPrint']) + "\n" +
             '       print: ' + str(option['print']) + "\n" +
             '  resultFile: ' + str(option['resultFile']) + "\n" +
             '  sdpaResult: ' + str(option['sdpaResult']) + "\n" +
             '  numThreads: ' + str(option['numThreads']) + "\n" +
             '   frvMethod: ' + str(option['frvMethod']) + "\n" +
             '  convMethod: ' + str(option['convMethod']) + "\n" +
             'domainMethod: ' + str(option['domainMethod']) + "\n" +
             ' rangeMethod: ' + str(option['rangeMethod']) + "\n" +
             '         rho: ' + str(option['rho']) + "\n" +
             '   zeroPoint: ' + str(option['zeroPoint']) + "\n\n")

    return

def write_symcone(fp, K, J=None):
    """Write SymCone to file

    Args:
      fp: File pointer
      K, J: SymCone
    """
    if J != None:
        fp.write("K =\n" + str(K) + "J =\n" + str(J) + "\n")
    else:
        fp.write("K =\n" + str(K) + "\n")

    return

def write_info(fp, sdpapinfo, sdpainfo, timeinfo):
    """Write result info to file

    Args:
      fp: File pointer
      sdpapinfo, sdpainfo, timeinfo: Result information
    """
    fp.write("----- Result -----\n" +
             "    SDPA.phase = %s\n" % sdpainfo['phasevalue'] +
             "     iteration = %d\n" % sdpainfo['iteration'] +
             "     primalObj = %+10.16e\n" % sdpapinfo['primalObj'] +
             "       dualObj = %+10.16e\n" % sdpapinfo['dualObj'] +
             "    dualityGap = %+10.16e\n" % sdpapinfo['dualityGap'] +
             "   primalError = %+10.16e\n" % sdpapinfo['primalError'] +
             "     dualError = %+10.16e\n" % sdpapinfo['dualError'] +
             "   convertTime = %f\n" % timeinfo['convert'] +
             "     solveTime = %f\n" % timeinfo['sdpa'] +
             "retrievingTime = %f\n" % timeinfo['retrieve'] +
             "     totalTime = %f\n" % timeinfo['total'])

    return

def write_result(fp, x, y):
    """Write result x and y to file

    Args:
      fp: File pointer
      x, y: results
    """
    if not sparse.isspmatrix_csc(x):
        x = x.tocsc()
    if not sparse.isspmatrix_csc(y):
        y = y.tocsc()

    x.sort_indices()
    y.sort_indices()

    fp.write("\nx(index, value) =\n")
    writeStr = ''
    count = 0
    for (index, value) in zip(x.indices, x.data):
        writeStr += "(%d: %+8.3e), " % (index, value)
        if count == 4:
            writeStr += "\n"
            count = 0
        else:
            count += 1

    fp.write(writeStr)

    fp.write("\n\ny(index, value) =\n")
    writeStr = ''
    count = 0
    for (index, value) in zip(y.indices, y.data):
        writeStr += "(%d: %+8.3e), " % (index, value)
        if count == 4:
            writeStr += "\n"
            count = 0
        else:
            count += 1


    fp.write(writeStr)

    return
