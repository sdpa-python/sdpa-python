#!/usr/bin/env python
"""sdpap.py
SDPAP (SDPA Python Interface)

This is the module of sdpap.
SDPAP solves a Conic-form Linear Optimization Problem (CLP) of the form:
minimize c^T x
subject to Ax - b \in J,  x \in K
K = (K.f, K.l, K.q, K.s)
J = (J.f, J.l, J.q, J.s)

September 2010, Kenta KATO

December 2010: modified for scipy
"""

__all__ = ['solve']

from . import convert
from . import sdpaputils
from . import fileio
from .param import param
from .symcone import SymCone
from .sdpacall import sdpacall
from .spcolo import spcolo
from .fvelim import fvelim
from scipy import sparse
from numpy import matrix
import copy
import time

def solve(A, b, c, K, J, option={}):
    """Solve CLP by SDPA

    If J.l or J.q or J.s > 0, clp_toLMI() or clp_toEQ() is called before solve.

    Args:
      A, b, c: Scipy matrices to denote the CLP.
      K, J: Symcone object to denote the CLP.
      option: Parameters. If None, default parameters is used.

    Returns:
      A tuple (x, y, sdpapinfo, timeinfo, sdpainfo).
      x, y: Primal and Dual solutions
      sdpapinfo, timeinfo, sdpainfo: Result information
    """
    timeinfo = dict()
    timeinfo['total'] = time.time()

    if 'print' not in option:
        option['print'] = 'display'
    verbose = len(option['print']) != 0 and option['print'] != 'no'
    maybe_print = print if verbose else lambda *a, **k: None

    # --------------------------------------------------
    # Set parameter
    # --------------------------------------------------
    option = param(option)
    maybe_print('---------- SDPAP Start ----------')

    # Write to output file
    if option['resultFile']:
        fpout = open(option['resultFile'], 'w')
        fileio.write_version(fpout)
        fileio.write_parameter(fpout, option)

    # --------------------------------------------------
    # Check validity
    # --------------------------------------------------
    if not K.check_validity() or not J.check_validity():
        return None

    if not isinstance(b, matrix) and not sparse.issparse(b):
        raise TypeError('sdpap.solve(): b must be a matrix or a sparse matrix.')
    if not isinstance(c, matrix) and not sparse.issparse(c):
        raise TypeError('sdpap.solve(): c must be a matrix or a sparse matrix.')
    if not isinstance(A, matrix) and not sparse.issparse(A):
        raise TypeError('sdpap.solve(): A must be a matrix or a sparse matrix.')

    if not sparse.isspmatrix_csc(b):
        b = sparse.csc_matrix(b)
    if b.shape[1] != 1:
        b = (b.T).tocsc()

    if not sparse.isspmatrix_csc(c):
        c = sparse.csc_matrix(c)
    if c.shape[1] != 1:
        c = (c.T).tocsc()

    if not sparse.isspmatrix_csc(A):
        A = sparse.csc_matrix(A)

    size_row = max(b.shape)
    size_col = max(c.shape)
    mA, nA = A.shape

    totalSize_n = K.f + K.l + sum(K.q) + sum(z ** 2 for z in K.s)
    totalSize_m = J.f + J.l + sum(J.q) + sum(z ** 2 for z in J.s)
    if size_row != mA or size_col != nA:
        maybe_print("Size A[m = %d, n = %d], b[m = %d], c[n = %d] ::" %
              (mA, nA, size_row, size_col))
        maybe_print("nnz(A) = %d, nnz(c) = %d" % (A.nnz, c.nnz))
        raise ValueError('Inconsistent Size')
    if size_col != totalSize_n:
        maybe_print("Size A[m = %d, n = %d], b[m = %d], c[n = %d] ::" %
              (mA, nA, size_row, size_col))
        maybe_print("nnz(A) = %d, nnz(c) = %d" % (A.nnz, c.nnz))
        raise ValueError("Inconsistent Size c[n = %d], K[%d]"
                         % (size_col, totalSize_n))
    if size_row != totalSize_m:
        maybe_print("Size A[m = %d, n = %d], b[m = %d], c[n = %d] ::" %
              (mA, nA, size_row, size_col))
        maybe_print("nnz(A) = %d, nnz(c) = %d" % (A.nnz, c.nnz))
        raise ValueError("Inconsistent Size b[n = %d], J[%d]"
                         % (size_row, totalSize_m))

    if option['resultFile']:
        fpout.write("----- Input Problem -----\n")
        fileio.write_symcone(fpout, K, J)

    # --------------------------------------------------
    # Exploiting sparsity conversion
    # --------------------------------------------------
    timeinfo['conv_domain'] = time.time()

    # Convert domain space sparsity
    if len(K.s) > 0:
        if option['domainMethod'] == 'clique':
            maybe_print('Applying the d-space conversion method '
                  'using clique trees...')
            dom_A, dom_b, dom_c, dom_K, dom_J, cliqueD = \
                   spcolo.dconv_cliquetree(A, b, c, K, J)
            ############################################################
            # Under construction
            ############################################################
            return
        elif option['domainMethod'] == 'basis':
            maybe_print('Applying the d-space conversion method '
                  'using basis representation...')
            dom_A, dom_b, dom_c, dom_K, dom_J, cliqueD = \
                   spcolo.dconv_basisrep(A, b, c, K, J)
        else:
            dom_A = copy.deepcopy(A)
            dom_b = copy.deepcopy(b)
            dom_c = copy.deepcopy(c)
            dom_K = copy.deepcopy(K)
            dom_J = copy.deepcopy(J)
    else:
        dom_A = copy.deepcopy(A)
        dom_b = copy.deepcopy(b)
        dom_c = copy.deepcopy(c)
        dom_K = copy.deepcopy(K)
        dom_J = copy.deepcopy(J)

    timeinfo['conv_domain'] = time.time() - timeinfo['conv_domain']

    if option['resultFile'] and option['domainMethod'] != 'none':
        fpout.write("----- Domain Space Sparsity Converted Problem-----\n")
        fileio.write_symcone(fpout, dom_K, dom_J)

    # Convert range space sparsity
    timeinfo['conv_range'] = time.time()
    if len(dom_J.s) > 0:
        if option['rangeMethod'] == 'clique':
            maybe_print('Applying the r-space conversion method '
                  'using clique trees...')
            ran_A, ran_b, ran_c, ran_K, ran_J, cliqueR = \
                    spcolo.rconv_cliquetree(dom_A, dom_b, dom_c, dom_K, dom_J)
            ############################################################
            # Under construction
            ############################################################
        elif option['rangeMethod'] == 'decomp':
            maybe_print('Applying the r-space conversion method '
                  'using matrix decomposition...')
            ran_A, ran_b, ran_c, ran_K, ran_J, cliqueR = \
                    spcolo.rconv_matdecomp(dom_A, dom_b, dom_c, dom_K, dom_J)
        else:
            ran_A, ran_b, ran_c = dom_A, dom_b, dom_c
            ran_K = copy.deepcopy(dom_K)
            ran_J = copy.deepcopy(dom_J)
    else:
        ran_A, ran_b, ran_c = dom_A, dom_b, dom_c
        ran_K = copy.deepcopy(dom_K)
        ran_J = copy.deepcopy(dom_J)

    timeinfo['conv_range'] = time.time() - timeinfo['conv_range']

    if option['resultFile'] and option['rangeMethod'] != 'none':
        fpout.write("----- Range Space Sparsity Converted Problem-----\n")
        fileio.write_symcone(fpout, ran_K, ran_J)

    # --------------------------------------------------
    # Convert to SeDuMi standard form
    # --------------------------------------------------
    timeinfo['conv_std'] = time.time()

    useConvert = False
    if ran_J.l > 0 or len(ran_J.q) > 0 or len(ran_J.s) > 0:
        useConvert = True
        if option['convMethod'] == 'LMI':
            maybe_print('Converting CLP format to LMI standard form...')
            A2, b2, c2, K2, J2, map_sdpIndex = \
                convert.clp_toLMI(ran_A, ran_b, ran_c, ran_K, ran_J)
        elif option['convMethod'] == 'EQ':
            maybe_print('Converting CLP format to EQ standard form.')
            ##################################################
            # This method is under construction
            ##################################################
            A2, b2, c2, K2, J2 = \
                convert.clp_toEQ(ran_A, ran_b, ran_c, ran_K, ran_J)
        else:
            raise ValueError("convMethod must be 'LMI' or 'EQ'")
    else:
        A2, b2, c2 = ran_A, ran_b, ran_c
        K2 = copy.deepcopy(ran_K)
        J2 = copy.deepcopy(ran_J)

    timeinfo['conv_std'] = time.time() - timeinfo['conv_std']

    if option['resultFile'] and \
           (ran_J.l > 0 or len(ran_J.q) > 0 or len(ran_J.s) > 0):
        fpout.write("----- SeDuMi format Converted Problem-----\n")
        fileio.write_symcone(fpout, K2)

    # --------------------------------------------------
    # Eliminate free variables
    # --------------------------------------------------
    timeinfo['conv_fv'] = time.time()

    if K2.f > 0:
        if option['frvMethod'] == 'split':
            maybe_print('Eliminating free variables with split method...')
            A3, b3, c3, K3 = fvelim.split(A2, b2, c2, K2, option['rho'])
        elif option['frvMethod'] == 'elimination':
            maybe_print('Eliminationg free variables with elimination method...')
            (A3, b3, c3, K3,
             LiP, U, Q, LPA_B, LPb_B, cfQU, gamma, rank_Af) = \
                fvelim.eliminate(A2, b2, c2, K2,
                                 option['rho'], option['zeroPoint'])
        else:
            raise ValueError("frvMethod must be 'split' or 'elimination'")
    else:
        A3, b3, c3, K3 = A2, b2, c2, K2

    timeinfo['conv_fv'] = time.time() - timeinfo['conv_fv']

    if option['resultFile'] and K2.f > 0:
        fpout.write("----- Free Variables Eliminated Problem -----\n")
        fileio.write_symcone(fpout, K3)

    # --------------------------------------------------
    # Solve by SDPA
    # --------------------------------------------------
    timeinfo['sdpa'] = time.time()
    x3, y3, s3, sdpainfo = sdpacall.solve_sdpa(A3, b3, c3, K3, option)
    timeinfo['sdpa'] = time.time() - timeinfo['sdpa']


    # --------------------------------------------------
    # Get Result
    # --------------------------------------------------
    maybe_print('Start: getCLPresult')

    # Retrieve result of fvelim
    timeinfo['ret_fv'] = time.time();
    if K2.f > 0:
        if option['frvMethod'] == 'split':
            maybe_print('Retrieving result with split method...')
            x2, y2, s2 = fvelim.result_split(x3, y3, s3, K2)
        elif option['frvMethod'] == 'elimination':
            maybe_print('Retrieving result with elimination method...')
            x2, y2, s2 = fvelim.result_elimination(x3, y3, s3, K2,
                                               LiP, U, Q,
                                               LPA_B, LPb_B, cfQU, rank_Af)
        else:
            raise ValueError("frvMethod must be 'split' or 'elimination'")
    else:
        x2, y2, s2 = x3, y3, s3

    timeinfo['ret_fv'] = time.time() - timeinfo['ret_fv']

    # Retrieve result from LMI or EQ
    timeinfo['ret_std'] = time.time()
    if useConvert:
        if option['convMethod'] == 'LMI':
            maybe_print('Retrieving result from LMI standard form...')
            x, y = convert.result_fromLMI(x2, y2, ran_K, ran_J, map_sdpIndex)
            tmp = -sdpainfo['primalObj']
            sdpainfo['primalObj'] = -sdpainfo['dualObj']
            sdpainfo['dualObj'] = tmp
        elif option['convMethod'] == 'EQ':
            maybe_print('Retrieving result from EQ standard form...')
            ##################################################
            # This method is under construction
            ##################################################
            x, y = result_fromEQ(x2, y2, ran_K, ran_J)
        else:
            raise ValueError("Something wrong about option['convMethod']")
    else:
        x, y = x2, y2

    timeinfo['ret_std'] = time.time() - timeinfo['ret_std']

    # Retrieve an optiomal solution from range space sparsity converted problem
    timeinfo['ret_range'] = time.time()
    if option['rangeMethod'] != 'none' and len(J.s) > 0:
        if option['rangeMethod'] == 'clique':
            maybe_print('Retrieving result with r-space conversion method '
                  'using clique trees...')
            ############################################################
            # Under construction
            ############################################################
            x, y = spcolo.rconv_cliqueresult(x, y, dom_K, dom_J, ran_K, cliqueR)
        elif option['rangeMethod'] == 'decomp':
            maybe_print('Retrieving result with r-space conversion method '
                  'using matrix decomposition...')
            x, y = spcolo.rconv_decompresult(x, y, dom_K, dom_J, ran_J, cliqueR)

    timeinfo['ret_range'] = time.time() - timeinfo['ret_range']

    # Retrieve an optiomal solution from domain space sparsity converted problem
    timeinfo['ret_domain'] = time.time()
    if option['domainMethod'] != 'none' and len(K.s) > 0:
        if option['domainMethod'] == 'clique':
            maybe_print('Retrieving result with d-space conversion method '
                  'using clique trees...')
            ############################################################
            # Under construction
            ############################################################
            x, y = spcolo.dconv_cliqueresult(x, y, K, J, dom_J, cliqueD)
        elif option['domainMethod'] == 'basis':
            maybe_print('Retrieving result with d-space conversion method '
                  'using basis representation...')
            x, y = spcolo.dconv_basisresult(x, y, K, J, dom_K, cliqueD)

    timeinfo['ret_domain'] = time.time() - timeinfo['ret_domain']
    timeinfo['total'] = time.time() - timeinfo['total']

    # --------------------------------------------------
    # Make dictionary 'info'
    # --------------------------------------------------
    maybe_print('Making result infomation...')
    timeinfo['convert'] = (timeinfo['conv_domain'] + timeinfo['conv_range'] +
                           timeinfo['conv_std'] + timeinfo['conv_fv'])
    timeinfo['retrieve'] = (timeinfo['ret_fv'] + timeinfo['ret_std'] +
                            timeinfo['ret_range'] + timeinfo['ret_domain'])

    if ran_K.f > 0 and option['frvMethod'] == 'elimination':
        sdpainfo['primalObj'] += gamma
        sdpainfo['dualObj'] += gamma

    sdpapinfo = dict()
    sdpapinfo['primalObj'] = (c.T * x)[0,0]
    sdpapinfo['dualObj'] = (b.T * y)[0,0]
    sdpapinfo['dualityGap'] = sdpaputils.get_dualitygap(x, y, b, c)
    sdpapinfo['primalError'] = sdpaputils.get_primalerror(x, A, b, J)
    sdpapinfo['dualError'] = sdpaputils.get_dualerror(y, A, c, K)

    # --------------------------------------------------
    # Print result
    # --------------------------------------------------
    if option['resultFile']:
        fileio.write_info(fpout, sdpapinfo, sdpainfo, timeinfo)
        fileio.write_result(fpout, x, y)
        fpout.close()

    maybe_print('========================================')
    maybe_print(' SDPAP: Result')
    maybe_print('========================================')
    maybe_print("    SDPA.phase = %s" % sdpainfo['phasevalue'])
    maybe_print("     iteration = %d" % sdpainfo['iteration'])
    maybe_print("    convMethod = %s" % option['convMethod'])
    maybe_print("     frvMethod = %s" % option['frvMethod'])
    maybe_print("  domainMethod = %s" % option['domainMethod'])
    maybe_print("   rangeMethod = %s" % option['rangeMethod'])
    maybe_print("     primalObj = %+10.16e" % sdpapinfo['primalObj'])
    maybe_print("       dualObj = %+10.16e" % sdpapinfo['dualObj'])
    maybe_print("    dualityGap = %+10.16e" % sdpapinfo['dualityGap'])
    maybe_print("   primalError = %+10.16e" % sdpapinfo['primalError'])
    maybe_print("     dualError = %+10.16e" % sdpapinfo['dualError'])
    maybe_print("   convertTime = %f" % timeinfo['convert'])
    maybe_print("     solveTime = %f" % timeinfo['sdpa'])
    maybe_print("retrievingTime = %f" % timeinfo['retrieve'])
    maybe_print("     totalTime = %f" % timeinfo['total'])
    maybe_print('---------- SDPAP End ----------')

    return x, y, sdpapinfo, timeinfo, sdpainfo
