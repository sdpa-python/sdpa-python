/*
Basic routines required by SparseCoLO subprograms
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

October 2010: Originally written by Kenta Kato
March 2022: Modified for Python 3 headers by Usama Muneeb
*/

#include <Python.h>
#define NO_ANSI99_COMPLEX
#include <iostream>
#include <cstdio>
#include <math.h>
#include "spcoloextmodule.h"
#include "spcolo_SparseMatrix.h"

// just for compatibility
#define mwSize Py_ssize_t
#define mwIndex Py_ssize_t

#define lengthOfString 10240
#define MX_DEBUG 0


/* ============================================================
   Get the permutation ordering by minimum degree ordering algorithm
   ============================================================ */
static char doc_ordering_mmd[] =
    "ordering = spcoloext.ordering_mmd(asp)";

static PyObject* ordering_mmd(PyObject* self, PyObject* args)
{
    PyObject* asp_ptr = NULL;

    if (!PyArg_ParseTuple(args, "O", &asp_ptr)) {
        Py_RETURN_NONE;
    }

#if MX_DEBUG
    rMessage("");
#endif

    //PyObject* module = PyImport_ImportModule("matdata.MatData");

    PyObject* tmpObj = NULL;

    /* Get size_n */
    tmpObj = PyObject_GetAttrString(asp_ptr, "size_col");
    int size_n = PyLong_AsLong(tmpObj);

#if MX_DEBUG
    printf("size_n = %d\n", size_n);
    rMessage("");
#endif

    /* Get num of nonzero */
    tmpObj = PyObject_GetAttrString(asp_ptr, "values");
    int num_nnz = PyList_Size(tmpObj);

#if MX_DEBUG
    rMessage("");
#endif

    /* Get values */
    double* asp_values = new double[num_nnz];
    for (int i = 0; i < num_nnz; i++) {
        asp_values[i] = PyFloat_AsDouble(PyList_GetItem(tmpObj, i));
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* Get rowind */
    tmpObj = PyObject_GetAttrString(asp_ptr, "rowind");
    int* asp_rowind = new int[num_nnz];
    for (int i = 0; i < num_nnz; i++) {
        asp_rowind[i] = PyLong_AsLong(PyList_GetItem(tmpObj, i));
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* Get colptr */
    tmpObj = PyObject_GetAttrString(asp_ptr, "colptr");
    int * asp_colptr = new int[size_n + 1];
    for (int i = 0; i <= size_n; i++) {
        asp_colptr[i] = PyLong_AsLong(PyList_GetItem(tmpObj, i));
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* Return object */
    PyObject* retObj = NULL;

    /* --------------------------------------------------
       Get ordering via SPOOLES
       -------------------------------------------------- */
    int* ordering = spcolo_ordering_mmd(asp_rowind, asp_colptr, size_n);

    delete [] asp_values;
    delete [] asp_rowind;
    delete [] asp_colptr;

    if (ordering == NULL) {
        /* Not need to decompose matrix */
        retObj = PyList_New(0);
        return retObj;
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       Get ordering list and copy to PyObject
       -------------------------------------------------- */
    retObj = PyList_New(size_n);
    for (int i = 0; i < size_n; i++) {
        PyList_SetItem(retObj, i, Py_BuildValue("i", ordering[i]));
    }

    return Py_BuildValue("O", retObj);
}



/* ============================================================
   Cholesky factorization
   ============================================================ */
static char doc_cholesky[] =
    "R = spcoloext.cholesky(A)";

static PyObject* cholesky(PyObject* self, PyObject* args)
{
    PyObject* A_ptr = NULL;

    if (!PyArg_ParseTuple(args, "O", &A_ptr)) {
        printf("scloext.cholesky(): Failed to partse tuple.\n");
        Py_RETURN_NONE;
    }

    PyObject* tmpObj = NULL;

    /* Get size_n */
    tmpObj = PyObject_GetAttrString(A_ptr, "size_col");
    int size_n = PyLong_AsLong(tmpObj);

    /* Get num of nonzero */
    tmpObj = PyObject_GetAttrString(A_ptr, "values");
    int num_nnz = PyList_Size(tmpObj);

    /* Get values */
    double* A_values = new double[num_nnz];
    for (int i = 0; i < num_nnz; i++) {
        A_values[i] = PyFloat_AsDouble(PyList_GetItem(tmpObj, i));
    }

    /* Get rowind */
    tmpObj = PyObject_GetAttrString(A_ptr, "rowind");
    int* A_rowind = new int[num_nnz];
    for (int i = 0; i < num_nnz; i++) {
        A_rowind[i] = PyLong_AsLong(PyList_GetItem(tmpObj, i));
    }

    /* Get colptr */
    tmpObj = PyObject_GetAttrString(A_ptr, "colptr");
    int * A_colptr = new int[size_n + 1];
    for (int i = 0; i <= size_n; i++) {
        A_colptr[i] = PyLong_AsLong(PyList_GetItem(tmpObj, i));
    }

    SparseMatrix* A = new SparseMatrix(size_n, size_n, num_nnz);

    unsigned int start_ptr, end_ptr;
    for (int col = 0; col < size_n; col++) {
        start_ptr = A_colptr[col];
        end_ptr = A_colptr[col + 1];
        for (int index = start_ptr; index < end_ptr; index++) {
            A->pushBack(A_rowind[index], col, A_values[index]);
        }
    }

    delete [] A_values;
    delete [] A_rowind;
    delete [] A_colptr;

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       Cholesky factorization
       -------------------------------------------------- */
    SparseMatrix* R = spcolo_cholesky(A);
    delete A;
    A = NULL;

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       Make return object
       -------------------------------------------------- */
    PyObject* valObj = PyList_New(R->getNumNonzero());
    PyObject* rowObj = PyList_New(R->getNumNonzero());
    PyObject* colObj = PyList_New(size_n + 1);
    if (!valObj || !rowObj || !colObj) {
        printf("Memory Over for Solution Space\n");
        Py_XDECREF(valObj);
        Py_XDECREF(rowObj);
        Py_XDECREF(colObj);
        return PyErr_NoMemory();
    }

    for (int i = 0; i < R->getNumNonzero(); i++) {
        PyList_SetItem(valObj, i, PyFloat_FromDouble(R->values[i]));
        PyList_SetItem(rowObj, i, PyLong_FromLong(R->rowind[i]));
    }

    for (int i = 0; i <= size_n; i++) {
        PyList_SetItem(colObj, i, PyLong_FromLong(R->colptr[i]));
    }

#if MX_DEBUG
    rMessage("");
#endif

    delete R;
    R = NULL;

#if MX_DEBUG
    rMessage("");
#endif


    return Py_BuildValue("OOO", valObj, rowObj, colObj);
}


/* ============================================================
   Initialize of this module
   ============================================================ */
PyDoc_STRVAR(spcoloext__doc__, "SPCOLOEXT: spcolo internal API.\n **** CAUTION **** \nDo NOT call directly. Call via spcoloext module.\n");

static PyObject* spcoloextmodule;

static PyMethodDef spcoloext_methods[] = {
    {"ordering_mmd", (PyCFunction)ordering_mmd, METH_VARARGS, doc_ordering_mmd},
    {"cholesky", (PyCFunction)cholesky, METH_VARARGS, doc_cholesky},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef spcoloext =
{
    PyModuleDef_HEAD_INIT,
    "spcoloext",
    NULL,
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    spcoloext_methods
};

PyMODINIT_FUNC PyInit_spcoloext(void)
{
    return PyModule_Create(&spcoloext);
}
