/* ************************************************************
   fvelimextmodule.cpp
   A extention module of fvelim with C or C++

   December 2010, Kenta KATO
   ************************************************************ */

#include <Python.h>
#define NO_ANSI99_COMPLEX
#include <iostream>
#include <cstdio>
#include <math.h>
#include "fvelim_SparseMatrix.h"
#include "fvelim_LUFactor.h"
#include "fvelim_util.h"

// just for compatibility
#define mwSize Py_ssize_t
#define mwIndex Py_ssize_t

#define lengthOfString 10240
#define MX_DEBUG 0


/* ============================================================
   LU factorization
   ============================================================ */
static char doc_lu[] =
    "LiP, U, P, Q, rank_Af = fvelimext.lu(Af, Alqs, rho, zeroPoint)";

static PyObject* lu(PyObject* self, PyObject* args)
{
    PyObject* Af_ptr = NULL;
    PyObject* Alqs_ptr = NULL;
    double rho;
    double zeroPoint;

    if (!PyArg_ParseTuple(args, "OOdd", &Af_ptr, &Alqs_ptr, &rho, &zeroPoint)) {
        Py_RETURN_NONE;
    }

    PyObject* tmpObj = NULL;
    int size_row, size_col, num_nnz, size_Kf;
    double* values = NULL;
    int* rowind = NULL;
    int* colptr = NULL;

    int num_ele;
    int ind_ptr = 0;

#if MX_DEBUG
    rMessage("");
#endif

    /* Make A_f */
    tmpObj = PyObject_GetAttrString(Af_ptr, "size_row");
    size_row = PyLong_AsLong(tmpObj);
    tmpObj = PyObject_GetAttrString(Af_ptr, "size_col");
    size_col = PyLong_AsLong(tmpObj);
    size_Kf = size_col;

    tmpObj = PyObject_GetAttrString(Af_ptr, "values");
    num_nnz = PyList_Size(tmpObj);
    NewArray(values, double, num_nnz);
    /*
    values = new double[num_nnz];
    */
    for (int i = 0; i < num_nnz; i++) {
        values[i] = PyFloat_AsDouble(PyList_GetItem(tmpObj, i));
    }

    tmpObj = PyObject_GetAttrString(Af_ptr, "rowind");
    NewArray(rowind, int, num_nnz);
    for (int i = 0; i < num_nnz; i++) {
        rowind[i] = PyLong_AsLong(PyList_GetItem(tmpObj, i));
    }

    tmpObj = PyObject_GetAttrString(Af_ptr, "colptr");
    NewArray(colptr, int, size_col + 1);
    for (int i = 0; i <= size_col; i++) {
        colptr[i] = PyLong_AsLong(PyList_GetItem(tmpObj, i));
    }

    SparseMatrix* A_f = new SparseMatrix(size_row, size_col, num_nnz);

    ind_ptr = 0;
    for (int col = 0; col < size_col; col++) {
        num_ele = colptr[col + 1] - colptr[col];
        for (int i = 0; i < num_ele; i++) {
            A_f->pushBack(rowind[ind_ptr], col, values[ind_ptr]);
            ind_ptr++;
        }
    }

    DeleteArray(values);
    DeleteArray(rowind);
    DeleteArray(colptr);

#if MX_DEBUG
    rMessage("");
#endif

    /* Make A_lqs */
    tmpObj = PyObject_GetAttrString(Alqs_ptr, "size_col");
    size_col = PyLong_AsLong(tmpObj);

    tmpObj = PyObject_GetAttrString(Alqs_ptr, "values");
    num_nnz = PyList_Size(tmpObj);
    NewArray(values, double, num_nnz);
    for (int i = 0; i < num_nnz; i++) {
        values[i] = PyFloat_AsDouble(PyList_GetItem(tmpObj, i));
    }

    tmpObj = PyObject_GetAttrString(Alqs_ptr, "rowind");
    NewArray(rowind, int, num_nnz);
    for (int i = 0; i < num_nnz; i++) {
        rowind[i] = PyLong_AsLong(PyList_GetItem(tmpObj, i));
    }

    tmpObj = PyObject_GetAttrString(Alqs_ptr, "colptr");
    NewArray(colptr, int, size_col + 1);
    for (int i = 0; i <= size_col; i++) {
        colptr[i] = PyLong_AsLong(PyList_GetItem(tmpObj, i));
    }

    SparseMatrix* A_lqs = new SparseMatrix(size_row, size_col, num_nnz);

    ind_ptr = 0;
    for (int col = 0; col < size_col; col++) {
        num_ele = colptr[col + 1] - colptr[col];
        for (int i = 0; i < num_ele; i++) {
            A_lqs->pushBack(rowind[ind_ptr], col, values[ind_ptr]);
            ind_ptr++;
        }
    }

    DeleteArray(values);
    DeleteArray(rowind);
    DeleteArray(colptr);

#if MX_DEBUG
    rMessage("");
#endif

    /* Make instance LUFactor */
    LUFactor LU(A_f, A_lqs, rho, zeroPoint);
    delete A_f;
    A_f = NULL;
    delete A_lqs;
    A_lqs = NULL;

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       LU factorization
       -------------------------------------------------- */
    LU.decompose();

#if MX_DEBUG
    rMessage("");
#endif

    /* --------------------------------------------------
       Make return object
       -------------------------------------------------- */
    SparseMatrix* LiP = LU.getLInvP();
    SparseMatrix* U = LU.getU();
    int* Q = LU.getQ();
    int rank_Af = LU.getRank();
    SparseMatrix* V = NULL;
    if (rank_Af < size_Kf) {
        V = LU.getV();
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* LiP */
    PyObject* LiP_val = PyList_New(LiP->getNumNonzero());
    PyObject* LiP_row = PyList_New(LiP->getNumNonzero());
    PyObject* LiP_col = PyList_New(size_row + 1);
    PyObject* LiPObj = PyTuple_New(3);
    if (!LiP_val || !LiP_row || !LiP_col || !LiPObj) {
        printf("Memory Over for Solution Space\n");
        Py_XDECREF(LiP_val);
        Py_XDECREF(LiP_row);
        Py_XDECREF(LiP_col);
        Py_XDECREF(LiPObj);
        return PyErr_NoMemory();
    }

    for (int i = 0; i < LiP->getNumNonzero(); i++) {
        PyList_SetItem(LiP_val, i, PyFloat_FromDouble(LiP->values[i]));
        PyList_SetItem(LiP_row, i, PyLong_FromLong(LiP->rowind[i]));
    }

    for (int i = 0; i <= size_row; i++) {
        PyList_SetItem(LiP_col, i, PyLong_FromLong(LiP->colptr[i]));
    }

    PyTuple_SetItem(LiPObj, 0, LiP_val);
    PyTuple_SetItem(LiPObj, 1, LiP_row);
    PyTuple_SetItem(LiPObj, 2, LiP_col);

#if MX_DEBUG
    rMessage("");
#endif

    /* U */
    PyObject* U_val = PyList_New(U->getNumNonzero());
    PyObject* U_row = PyList_New(U->getNumNonzero());
    PyObject* U_col = PyList_New(rank_Af + 1);
    PyObject* UObj = PyTuple_New(3);
    if (!U_val || !U_row || !U_col || !UObj) {
        printf("Memory Over for Solution Space\n");
        Py_XDECREF(U_val);
        Py_XDECREF(U_row);
        Py_XDECREF(U_col);
        Py_XDECREF(UObj);
        return PyErr_NoMemory();
    }

    for (int i = 0; i < U->getNumNonzero(); i++) {
        PyList_SetItem(U_val, i, PyFloat_FromDouble(U->values[i]));
        PyList_SetItem(U_row, i, PyLong_FromLong(U->rowind[i]));
    }

    for (int i = 0; i <= rank_Af; i++) {
        PyList_SetItem(U_col, i, PyLong_FromLong(U->colptr[i]));
    }

    PyTuple_SetItem(UObj, 0, U_val);
    PyTuple_SetItem(UObj, 1, U_row);
    PyTuple_SetItem(UObj, 2, U_col);

#if MX_DEBUG
    rMessage("");
#endif

    /* V */
    PyObject* VObj;
    if (rank_Af < size_Kf) {
        PyObject* V_val = PyList_New(V->getNumNonzero());
        PyObject* V_row = PyList_New(V->getNumNonzero());
        PyObject* V_col = PyList_New(size_Kf - rank_Af + 1);
        VObj = PyTuple_New(3);
        if (!V_val || !V_row || !V_col || !VObj) {
            printf("Memory Over for Solution Space\n");
            Py_XDECREF(V_val);
            Py_XDECREF(V_row);
            Py_XDECREF(V_col);
            Py_XDECREF(VObj);
            return PyErr_NoMemory();
        }

        for (int i = 0; i < V->getNumNonzero(); i++) {
            PyList_SetItem(V_val, i, PyFloat_FromDouble(V->values[i]));
            PyList_SetItem(V_row, i, PyLong_FromLong(V->rowind[i]));
        }

        int size = size_Kf - rank_Af;
        for (int i = 0; i <= size; i++) {
            PyList_SetItem(V_col, i, PyLong_FromLong(V->colptr[i]));
        }

        PyTuple_SetItem(VObj, 0, V_val);
        PyTuple_SetItem(VObj, 1, V_row);
        PyTuple_SetItem(VObj, 2, V_col);
    } else {
        VObj = Py_None;
    }

#if MX_DEBUG
    rMessage("");
#endif

    /* Q */
    PyObject* QObj = PyList_New(rank_Af);
    if (!QObj) {
        printf("Memory Over for Solution Space\n");
        Py_XDECREF(QObj);
        return PyErr_NoMemory();
    }

    for (int i = 0; i < rank_Af; i++) {
        PyList_SetItem(QObj, i, PyLong_FromLong(Q[i]));
    }

#if MX_DEBUG
    rMessage("");
#endif

    return Py_BuildValue("OOOOi", LiPObj, UObj, VObj, QObj, rank_Af);
}


/* ============================================================
   Initialize of this module
   ============================================================ */
PyDoc_STRVAR(fvelimext__doc__, "FVELIMEXT: fvelim internal API.\n **** CAUTION **** \nDo NOT call directly. Call via fvelimext module.\n");

static PyObject* fvelimextmodule;

static PyMethodDef fvelimext_methods[] = {
    {"lu", (PyCFunction)lu, METH_VARARGS, doc_lu},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef fvelimext =
{
    PyModuleDef_HEAD_INIT,
    "fvelimext",
    NULL,
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    fvelimext_methods
};

PyMODINIT_FUNC PyInit_fvelimext(void)
{
    return PyModule_Create(&fvelimext);
}
