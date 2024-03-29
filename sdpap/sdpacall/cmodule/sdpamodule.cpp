/*
C Extension to call SDPA C API
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
March 2022: Modified for Python 3 headers by Usama Muneeb
*/

#include <Python.h>
#include <iostream>
#include <cstdio>
#include <algorithm>
using namespace std;
#include <sdpa_call.h>
using namespace sdpa;
using sdpa::Time;

// just for compatibility
#define mwSize Py_ssize_t
#define mwIndex Py_ssize_t

#define lengthOfString 10240
#define MX_DEBUG 0

#if USEGMP
#include <gmpxx.h>
#endif
/* ============================================================
   Message
   ============================================================ */
#define rMessage(message)                       \
    {cout << message << " :: line " << __LINE__ \
          << " in " << __FILE__ << endl; }

#define rError(message)                         \
    {cout << message << " :: line " << __LINE__ \
          << " in " << __FILE__ << endl;        \
        exit(false);}

/* ============================================================
   Allocate array
   ============================================================ */
#if 1
#define NewArray(val,type,number) \
  {val = NULL; \
    try{ val = new type[number]; } \
    catch(bad_alloc){ \
        rMessage("Memory Exhausted (bad_alloc)"); abort(); } \
    catch(...){ \
        rMessage("Fatal Error (related memory allocation"); abort(); } \
  }
#else
#define NewArray(val,type,number) \
  {rMessage("New Invoked"); \
   val = NULL; val = new type[number]; \
   if  (val==NULL) {rError("Over Memory");} \
  }
#endif

#define DeleteArray(val) \
  { if  (val!=NULL) { \
      delete[] val; \
      val = NULL; \
    } \
  }


static PyObject* get_backend_info(PyObject *self, PyObject *args)
{
    return Py_BuildValue("{s:i}", "gmp", USEGMP);
}

/* ============================================================
   sdpamodule.sedumiwrap(At, b, c, dictK, option)
   ============================================================ */
static char doc_sedumiwrap[] =
    "[x,y,s,info] = sdpa.sedumiwrap(At,b,c,dictK,option)";

static PyObject* sedumiwrap(PyObject* self, PyObject* args, PyObject* kwrds)
{
    char* kwlist[] = {(char*)"At", (char*)"b", (char*)"c", (char*)"dictK", (char*)"option", NULL};

    PyObject* At_ptr = NULL;
    PyObject* b_ptr = NULL;
    PyObject* c_ptr = NULL;
    PyDictObject* dictK_ptr = NULL;
    PyDictObject* option_ptr = NULL;

    if(!PyArg_ParseTupleAndKeywords(args, kwrds, "OOOOO", kwlist,
                                    &At_ptr, &b_ptr, &c_ptr, &dictK_ptr, &option_ptr)){
        Py_RETURN_NONE;
    }

    time_t ltime;
    time(&ltime);
    char string_time[1024];
    strcpy(string_time, ctime(&ltime));
    string_time[strlen(string_time) - 1] = '\0';

    PyObject* tmpObj;

#if USEGMP
    /* Get numerical precision from solver options */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "mpfPrecision");
    int mpfPrecision = (int)PyLong_AsLong(tmpObj);
    mpf_set_default_prec(mpfPrecision);
#endif
    SDPA sdpa;

    int maxIteration = 0;
    int nSymmChk = 0;
    int nDimacs = 0;

    /* strings for phase value */
    const char* szPhase[] = {
        "noINFO", "pFEAS", "dFEAS", "pdFEAS", "pdINF",
        "pFEAS_dINF", "pINF_dFEAS", "pdOPT", "pUNBD", "dUNBD"};

    /* output file */
    char* outfile = NULL;
    FILE* fp = NULL;
    FILE* fpResult = NULL;
    int nOutfile = 0;

    mwSize mDim;
    mwSize nBlock;

    /* temporary variables */
    mwIndex k;
    int size;
#if USEGMP
    mpf_class* tmp_ptr = NULL;
#else
    double* tmp_ptr = NULL;
#endif
    char* tmpPrint = NULL;

    TimeStart(SDPA_START);
    TimeStart(SDPA_CONVERT_START);

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       Set SDPA parameters by OPTIONS
       -------------------------------------------------- */
    /* Max Iteration */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "maxIteration");
    maxIteration = (int)PyLong_AsLong(tmpObj);
    sdpa.setParameterMaxIteration(maxIteration);

    /* epsilonStar */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "epsilonStar");
    sdpa.setParameterEpsilonStar(PyFloat_AsDouble(tmpObj));

    /* lambdaStar */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "lambdaStar");
    sdpa.setParameterLambdaStar(PyFloat_AsDouble(tmpObj));

    /* omegaStar */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "omegaStar");
    sdpa.setParameterOmegaStar(PyFloat_AsDouble(tmpObj));

    /* lowerBound */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "lowerBound");
    sdpa.setParameterLowerBound(PyFloat_AsDouble(tmpObj));

    /* upperBound */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "upperBound");
    sdpa.setParameterUpperBound(PyFloat_AsDouble(tmpObj));

    /* betaStar */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "betaStar");
    sdpa.setParameterBetaStar(PyFloat_AsDouble(tmpObj));

    /* betaBar */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "betaBar");
    sdpa.setParameterBetaBar(PyFloat_AsDouble(tmpObj));

    /* gammaStar */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "gammaStar");
    sdpa.setParameterGammaStar(PyFloat_AsDouble(tmpObj));

    /* epsilonDash */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "epsilonDash");
    sdpa.setParameterEpsilonDash(PyFloat_AsDouble(tmpObj));

    /* isSymmetric */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "isSymmetric");
    nSymmChk = (int)PyLong_AsLong(tmpObj);

    /* isDimacs */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "isDimacs");
    nDimacs = (int)PyLong_AsLong(tmpObj);

    /* yPrint */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "yPrint");
    sdpa.setParameterPrintXVec((char*)PyUnicode_AsUTF8(tmpObj));

    /* sPrint */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "sPrint");
    sdpa.setParameterPrintXMat((char*)PyUnicode_AsUTF8(tmpObj));

    /* xPrint */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "xPrint");
    sdpa.setParameterPrintYMat((char*)PyUnicode_AsUTF8(tmpObj));

    /* infPrint */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "infPrint");
    sdpa.setParameterPrintInformation((char*)PyUnicode_AsUTF8(tmpObj));

    /* print */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "print");
    const char* outfile2 = PyUnicode_AsUTF8(tmpObj);
    /* default setting is displaying information to stdout */
    fp = stdout;
    if (strlen(outfile2) == 0) {
        fp = NULL;
    } else {
        if (strncmp("display", outfile2, strlen(outfile2)) == 0) {
            fp = stdout;
        } else if (strncmp("no", outfile2, strlen(outfile2)) == 0) {
            fp = NULL;
        } else {
            fp = fopen(outfile2, "at");
            if (fp == NULL) {
                printf("Failed to open %s\n", outfile2);
                fp = stdout;
            } else {
                nOutfile = 1;
            }
        }
    }
    sdpa.setDisplay(fp);

    /* resultFile */
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "sdpaResult");
    const char* outfile3 = PyUnicode_AsUTF8(tmpObj);
    fpResult = NULL;
    if (strlen(outfile3) > 0) {
        if(strncmp("no", outfile3, strlen(outfile3)) == 0) {
            // printf("resultFile is NULL\n");
        } else {
            fpResult = fopen(outfile3, "w");
            if (fpResult == NULL) {
                printf("Failed to open %s\n", outfile3);
                printf("Skip the detail file\n");
            }
        }
    }
    sdpa.setResultFile(fpResult);

    if (fp) {
#if USEGMP
        fprintf( fp,"SDPA Multiprecision start at [%s]\n",string_time );
#else
        fprintf( fp,"SDPA start at [%s]\n",string_time );
#endif
    }
    if (fpResult) {
#if USEGMP
        fprintf( fpResult,"SDPA Multiprecision start at [%s]\n",string_time );
#else
        fprintf( fpResult,"SDPA start at [%s]\n",string_time );
#endif
    }

    /* numThreads */
#if USEGMP
#else
    tmpObj = PyDict_GetItemString((PyObject*)option_ptr, "numThreads");
    sdpa.setNumThreads((int)PyLong_AsLong(tmpObj));
#endif

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       initialize SDPA class members
       -------------------------------------------------- */
    /* mDim -!- ATTENTION: Matrix A is transposed -!- */
    tmpObj = PyObject_GetAttrString(At_ptr, "size_col");
    mDim = PyLong_AsLong(tmpObj);
    tmpObj = PyObject_GetAttrString(At_ptr, "size_row");
    int size_n = PyLong_AsLong(tmpObj);

    PyObject* x_ptr = NULL;
    PyObject* y_ptr = NULL;
    PyObject* s_ptr = NULL;
    PyDictObject* info_ptr = NULL;

    x_ptr = PyList_New(size_n);
    y_ptr = PyList_New(mDim);
    s_ptr = PyList_New(size_n);
    info_ptr = (PyDictObject*)PyDict_New();
    if (x_ptr == NULL || y_ptr == NULL || s_ptr == NULL || info_ptr == NULL) {
        cout << "Memory Over for Solution Space" << endl;
        Py_XDECREF(x_ptr);
        Py_XDECREF(y_ptr);
        Py_XDECREF(s_ptr);
        Py_XDECREF(info_ptr);
        return PyErr_NoMemory();
    }

    sdpa.inputConstraintNumber(mDim);

#if MX_DEBUG
    rMessage("");
#endif

    /* nBlock */
    nBlock = 0;

    /* get K.f */
    size_t K_f = 0;
    int isK_f = 0;
    tmpObj = PyDict_GetItemString((PyObject*)dictK_ptr, "f");
    K_f = (int)PyLong_AsLong(tmpObj);
    if (K_f > 0) {
        rError("SDPA does not support K.f");
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* get K.l */
    size_t K_l = 0;
    int isK_l = 0;
    tmpObj = PyDict_GetItemString((PyObject*)dictK_ptr, "l");
    K_l = (int)PyLong_AsLong(tmpObj);
    if (K_l > 0) {
        isK_l = 1;
        nBlock++;
    }


#if MX_DEBUG
    rMessage("");
#endif

    /* get K.q */
    tmpObj = PyDict_GetItemString((PyObject*)dictK_ptr, "q");
    int* K_q = NULL;
    int* K_socpConeStart = NULL;
    int K_socpNoCones = 0;

    K_socpNoCones = (int)PyTuple_Size(tmpObj);
    if (K_socpNoCones > 0) {
        NewArray(K_q, int, K_socpNoCones);
        NewArray(K_socpConeStart, int, K_socpNoCones + 1);
        K_socpConeStart[0] = 0;

        for (int block = 0; block < K_socpNoCones; block++) {
            PyObject* tmpInt = PyTuple_GetItem(tmpObj, block);
            K_q[block] = (int)PyLong_AsLong(tmpInt);
            K_socpConeStart[block + 1] = K_socpConeStart[block] + K_q[block];
            nBlock++;
        }
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* get K.s */
    tmpObj = PyDict_GetItemString((PyObject*)dictK_ptr, "s");
    int* K_s = NULL;
    int* K_sdpConeStart = NULL;
    int K_sdpNoCones = 0;

    K_sdpNoCones = (int)PyTuple_Size(tmpObj);
    if (K_sdpNoCones > 0) {
        NewArray(K_s, int, K_sdpNoCones);
        NewArray(K_sdpConeStart, int, K_sdpNoCones + 1);
        K_sdpConeStart[0] = 0;

        for (int block = 0; block < K_sdpNoCones; block++) {
            PyObject* tmpInt = PyTuple_GetItem(tmpObj, block);
            K_s[block] = (int)PyLong_AsLong(tmpInt);
            K_sdpConeStart[block + 1] = K_sdpConeStart[block] + K_s[block] * K_s[block];
            nBlock++;
        }
    }

    sdpa.inputBlockNumber(nBlock);

    int size_Kq = 0;
    int size_Ks = 0;

#if MX_DEBUG
    rMessage("");
#endif

    if (isK_l > 0) {
        sdpa.inputBlockSize(1, K_l);
        sdpa.inputBlockType(1, SDPA::LP);
    }

#if MX_DEBUG
    rMessage("");
#endif

    for (int block = 0; block < K_socpNoCones; block++) {
        sdpa.inputBlockSize(block + isK_l + 1, K_q[block]);
        sdpa.inputBlockType(block + isK_l + 1, SDPA::SOCP);
        size_Kq += K_q[block];
    }

#if MX_DEBUG
    rMessage("");
#endif

    for (int block = 0; block < K_sdpNoCones; block++) {
        sdpa.inputBlockSize(block + isK_l + K_socpNoCones + 1, K_s[block]);
        sdpa.inputBlockType(block + isK_l + K_socpNoCones + 1, SDPA::SDP);
        size_Ks += K_s[block] * K_s[block];
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       Execute initializeUpperTriangleSpace()
       -------------------------------------------------- */
    sdpa.initializeUpperTriangleSpace();

#if MX_DEBUG
    rMessage("");
#endif

    int numNonzero = 0;
    int num_nnz;
    int row, col;
    double val;
    PyObject* valObj;
    PyObject* rowObj;
    PyObject* colObj;
    /* --------------------------------------------------
       Input b
       -------------------------------------------------- */
    valObj = PyObject_GetAttrString(b_ptr, "values");
    rowObj = PyObject_GetAttrString(b_ptr, "rowind");
    num_nnz = PyList_Size(valObj);
    for (int i = 0; i < num_nnz; i++) {
        val = PyFloat_AsDouble(PyList_GetItem(valObj, i));
        row = PyLong_AsLong(PyList_GetItem(rowObj, i));
        sdpa.inputCVec(row + 1, -val);
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       Input c
       -------------------------------------------------- */
    valObj = PyObject_GetAttrString(c_ptr, "values");
    rowObj = PyObject_GetAttrString(c_ptr, "rowind");
    num_nnz = PyList_Size(valObj);
    int currentSocpCone = 0;
    int currentSdpCone = 0;
    for (int index = 0; index < num_nnz; index++) {
        val = PyFloat_AsDouble(PyList_GetItem(valObj, index));
        row = PyLong_AsLong(PyList_GetItem(rowObj, index));
        if (row < K_l) {
            sdpa.inputElement(0, 1, row + 1, row + 1, -val);
        } else if (row - K_l < size_Kq) {
            row -= K_l;
            while (K_socpConeStart[currentSocpCone + 1] <= row) {
                currentSocpCone++;
            }
            row -= K_socpConeStart[currentSocpCone];
            sdpa.inputElement(0, isK_l + currentSocpCone + 1, row + 1, row + 1, -val);
        } else {
            row -= (K_l + size_Kq);
            while (K_sdpConeStart[currentSdpCone + 1] <= row) {
                currentSdpCone++;
            }
            row -= K_sdpConeStart[currentSdpCone];
            int i = row / K_s[currentSdpCone];
            int j = row % K_s[currentSdpCone];
            if (i <= j) {
                sdpa.inputElement(0, isK_l + K_socpNoCones + currentSdpCone + 1, i + 1, j + 1, -val);
            }
        }
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       Input A
       -!- ATTENTION: Matrix A is transposed -!-
       -------------------------------------------------- */
    valObj = PyObject_GetAttrString(At_ptr, "values");
    colObj = PyObject_GetAttrString(At_ptr, "rowind");
    rowObj = PyObject_GetAttrString(At_ptr, "colptr");

    int row_ptr[mDim + 1];
    for (int i = 0; i <= mDim; i++) {
        row_ptr[i] = PyLong_AsLong(PyList_GetItem(rowObj, i));
    }

    int start_idx, end_idx;
    for (int row = 0; row < mDim; row++) {
        num_nnz = row_ptr[row + 1] - row_ptr[row];
        start_idx = row_ptr[row];
        end_idx = row_ptr[row + 1];
        currentSocpCone = 0;
        currentSdpCone = 0;
        for (int index = start_idx; index < end_idx; index++) {
            val = PyFloat_AsDouble(PyList_GetItem(valObj, index));
            col = PyLong_AsLong(PyList_GetItem(colObj, index));
            if (col < K_l) {
                sdpa.inputElement(row + 1, 1, col + 1, col + 1, -val);
            } else if (col - K_l < size_Kq) {
                col -= K_l;
                while (K_socpConeStart[currentSocpCone + 1] <= col) {
                    currentSocpCone++;
                }
                col -= K_socpConeStart[currentSocpCone];
                sdpa.inputElement(row + 1, isK_l + currentSocpCone + 1, col + 1, col + 1, -val);
            } else {
                col -= (K_l + size_Kq);
                while (K_sdpConeStart[currentSdpCone + 1] <= col) {
                    currentSdpCone++;
                }
                col -= K_sdpConeStart[currentSdpCone];
                int i = col / K_s[currentSdpCone];
                int j = col % K_s[currentSdpCone];
                if (i <= j) {
                    sdpa.inputElement(row + 1, isK_l + K_socpNoCones + currentSdpCone + 1,
                                      i + 1, j + 1, -val);
                }
            }
        }
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       Check the consistence of F, c
       -------------------------------------------------- */
    if (nSymmChk) {
        sdpa.initializeUpperTriangle(true);
    } else {
        sdpa.initializeUpperTriangle(false);
    }

#if MX_DEBUG
    rMessage("");
    //sdpa.writeInputSparse((char*)"b.dat-s", (char*)"%e");
#endif

    /* --------------------------------------------------
       Solve by SDPA
       -------------------------------------------------- */
#if MX_DEBUG
    rMessage("");
#endif

    sdpa.initializeSolve();

#if MX_DEBUG
    rMessage("");
#endif

    if (fp) {
        fprintf(fp, "Converted to SDPA internal data / ");
        fprintf(fp, "Starting SDPA main loop\n");
    }
    TimeEnd(SDPA_CONVERT_END);
    TimeStart(SDPA_SOLVE_START);
    sdpa.solve();
    TimeEnd(SDPA_SOLVE_END);

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       Set output values to arguments
       -------------------------------------------------- */
    TimeStart(SDPA_RETRIEVE_START);
    if (fp) {
        fprintf(fp, "Converting optimal solution to CLP format\n");
    }

    double read_double;

    /* Optimal value for xVec */
    tmp_ptr = sdpa.getResultXVec();
    if (tmp_ptr != NULL) {
        for (k = 0; k < mDim; k++) {
#if USEGMP
            read_double = mpf_get_d(tmp_ptr[k].get_mpf_t());
#else
            read_double = tmp_ptr[k];
#endif
            PyList_SetItem(y_ptr, k, PyFloat_FromDouble(read_double));
        }
    }

    /* Optimal value for YMat */
    if (isK_l > 0) {
        size = sdpa.getBlockSize(1);
        tmp_ptr = sdpa.getResultYMat(1);
        for (int index = 0; index < size; index++) {
#if USEGMP
            read_double = mpf_get_d(tmp_ptr[index].get_mpf_t());
#else
            read_double = tmp_ptr[index];
#endif
            PyList_SetItem(x_ptr, index, PyFloat_FromDouble(read_double));
        }
    }

    for (int l = 0; l < K_socpNoCones; l++) {
        size = sdpa.getBlockSize(l + isK_l + 1);
        tmp_ptr = sdpa.getResultYMat(l + isK_l + 1);
        for (int index = 0; index < size; index++) {
#if USEGMP
            read_double = mpf_get_d(tmp_ptr[index].get_mpf_t());
#else
            read_double = tmp_ptr[index];
#endif
            PyList_SetItem(x_ptr, K_l + K_socpConeStart[l] + index,
                           PyFloat_FromDouble(read_double));
        }
    }

    for (int l = 0; l < K_sdpNoCones; l++) {
        size = sdpa.getBlockSize(l + isK_l + K_socpNoCones + 1);
        tmp_ptr = sdpa.getResultYMat(l + isK_l + K_socpNoCones + 1);
        for (int index = 0; index < size * size; index++) {
#if USEGMP
            read_double = mpf_get_d(tmp_ptr[index].get_mpf_t());
#else
            read_double = tmp_ptr[index];
#endif
            PyList_SetItem(x_ptr, K_l + size_Kq + K_sdpConeStart[l] + index,
                           PyFloat_FromDouble(read_double));
        }
    }

    /* Optimal value for XMat */
    if (isK_l > 0) {
        size = sdpa.getBlockSize(1);
        tmp_ptr = sdpa.getResultXMat(1);
        for (int index = 0; index < size; index++) {
#if USEGMP
            read_double = mpf_get_d(tmp_ptr[index].get_mpf_t());
#else
            read_double = tmp_ptr[index];
#endif
            PyList_SetItem(s_ptr, index, PyFloat_FromDouble(read_double));
        }
    }

    for (int l = 0; l < K_socpNoCones; l++) {
        size = sdpa.getBlockSize(l + isK_l + 1);
        tmp_ptr = sdpa.getResultXMat(l + isK_l + 1);
        for (int index = 0; index < size; index++) {
#if USEGMP
            read_double = mpf_get_d(tmp_ptr[index].get_mpf_t());
#else
            read_double = tmp_ptr[index];
#endif
            PyList_SetItem(s_ptr, K_l + K_socpConeStart[l] + index,
                           PyFloat_FromDouble(read_double));
        }
    }

    for (int l = 0; l < K_sdpNoCones; l++) {
        size = sdpa.getBlockSize(l + isK_l + K_socpNoCones + 1);
        tmp_ptr = sdpa.getResultXMat(l + isK_l + K_socpNoCones + 1);
        for (int index = 0; index < size * size; index++) {
#if USEGMP
            read_double = mpf_get_d(tmp_ptr[index].get_mpf_t());
#else
            read_double = tmp_ptr[index];
#endif
            PyList_SetItem(s_ptr, K_l + size_Kq + K_sdpConeStart[l] + index,
                           PyFloat_FromDouble(read_double));
        }
    }

    TimeEnd(SDPA_RETRIEVE_END);

#if USEGMP
#else
    /* Dimacs Error Information */
    if (nDimacs != 0) {
        printf("Computing Dimacs Error\n");
        double dimacs_error[7];
        sdpa.getDimacsError(dimacs_error);
        PyDict_SetItemString((PyObject*)info_ptr, "dimacs",
                             Py_BuildValue("(dddddd)",
                                           dimacs_error[1],
                                           dimacs_error[2],
                                           dimacs_error[3],
                                           dimacs_error[4],
                                           dimacs_error[5],
                                           dimacs_error[6]));
    }
#endif

    /* Phase information */
    PyDict_SetItemString((PyObject*)info_ptr, "phasevalue",
                         Py_BuildValue( "s", szPhase[sdpa.getPhaseValue()]));

    /* Iteration */
    PyDict_SetItemString((PyObject*)info_ptr, "iteration",
                         Py_BuildValue("i", sdpa.getIteration()));

    /* primalObj */
#if USEGMP
            read_double = mpf_get_d(sdpa.getDualObj().get_mpf_t());
#else
            read_double = sdpa.getDualObj();
#endif
    PyDict_SetItemString((PyObject*)info_ptr, "primalObj",
                         Py_BuildValue("d", -read_double));

    /* dualObj */
#if USEGMP
            read_double = mpf_get_d(sdpa.getPrimalObj().get_mpf_t());
#else
            read_double = sdpa.getPrimalObj();
#endif
    PyDict_SetItemString((PyObject*)info_ptr, "dualObj",
                         Py_BuildValue("d", -read_double));

    /* primalError */
#if USEGMP
            read_double = mpf_get_d(sdpa.getDualError().get_mpf_t());
#else
            read_double = sdpa.getDualError();
#endif
    PyDict_SetItemString((PyObject*)info_ptr, "primalError",
                         Py_BuildValue("d", read_double));

    /* dualError */
#if USEGMP
            read_double = mpf_get_d(sdpa.getPrimalError().get_mpf_t());
#else
            read_double = sdpa.getPrimalError();
#endif
    PyDict_SetItemString((PyObject*)info_ptr, "dualError",
                         Py_BuildValue("d", read_double));

    /* digits */
    PyDict_SetItemString((PyObject*)info_ptr, "digits",
                         Py_BuildValue("d", sdpa.getDigits()));

    /* dualityGap */
#if USEGMP
            read_double = mpf_get_d(sdpa.getDualityGap().get_mpf_t());
#else
            read_double = sdpa.getDualityGap();
#endif
    PyDict_SetItemString((PyObject*)info_ptr, "dualityGap",
                         Py_BuildValue("d", read_double));

    /* mu */
#if USEGMP
            read_double = mpf_get_d(sdpa.getMu().get_mpf_t());
#else
            read_double = sdpa.getMu();
#endif
    PyDict_SetItemString((PyObject*)info_ptr, "mu",
                         Py_BuildValue("d", read_double));

    /* solveTime */
    PyDict_SetItemString((PyObject*)info_ptr, "solveTime",
                         Py_BuildValue("d", TimeCal(SDPA_SOLVE_START, SDPA_SOLVE_END)));

    /* convertingTime */
    PyDict_SetItemString((PyObject*)info_ptr, "convertingTime",
                         Py_BuildValue("d", TimeCal(SDPA_CONVERT_START, SDPA_CONVERT_END)));

    /* retrivingTime */
    PyDict_SetItemString((PyObject*)info_ptr, "retrievingTime",
                         Py_BuildValue("d", TimeCal(SDPA_RETRIEVE_START, SDPA_RETRIEVE_END)));

    /* totalTime */
    TimeEnd(SDPA_END);
    PyDict_SetItemString((PyObject*)info_ptr, "sdpaTime",
                         Py_BuildValue("d", TimeCal(SDPA_START,SDPA_END)));

    time(&ltime);
    strcpy(string_time, ctime(&ltime));
    string_time[strlen(string_time) - 1] = '\0';

    if (fp) {
#if USEGMP
        fprintf(fp,"SDPA Multiprecision end at [%s]\n", string_time);
#else
        fprintf(fp,"SDPA end at [%s]\n", string_time);
#endif
    }
    if (fpResult) {
#if USEGMP
        fprintf(fpResult,"SDPA Multiprecision end at [%s]\n", string_time);
#else
        fprintf(fpResult,"SDPA end at [%s]\n", string_time);
#endif
    }

    /* close output file */
    if (nOutfile) {
        fclose(fp);
    }

#if MX_DEBUG
    rMessage("");
#endif

    /*** Free allocated memory ****/
    if (K_socpNoCones > 0) {
        DeleteArray(K_q);
        DeleteArray(K_socpConeStart);
    }

#if MX_DEBUG
    rMessage("");
#endif

    if (K_sdpNoCones > 0) {
        DeleteArray(K_s);
        DeleteArray(K_sdpConeStart);
    }

#if MX_DEBUG
    rMessage("");
#endif

/*
    This is automatically called by the destructor of the SDPA class.
    Calling twice results in problems with the multiprecision backend.
    sdpa.terminate();
*/

#if MX_DEBUG
    rMessage("");
#endif

    return Py_BuildValue("OOOO", x_ptr, y_ptr, s_ptr, info_ptr);
}


/* --------------------------------------------------
   INITIALIZER of this module
   -------------------------------------------------- */
PyDoc_STRVAR(sdpa__doc__, "SDPA: SDPAP internal API.\n **** CAUTION **** \nDo NOT call them directly. Call via SDPAP module.\n");

static PyObject* sdpamodule;

static PyMethodDef sdpa_methods[] = {
    {"get_backend_info", (PyCFunction)get_backend_info, METH_VARARGS, "Returns backend information."},
    {"sedumiwrap", (PyCFunction) sedumiwrap, METH_VARARGS | METH_KEYWORDS, doc_sedumiwrap},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef _sdpa =
{
    PyModuleDef_HEAD_INIT,
    "sdpa",
    NULL,
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    sdpa_methods
};

PyMODINIT_FUNC PyInit_sdpa(void)
{
    return PyModule_Create(&_sdpa);
}
