/* ************************************************************
   sclo_cholesky.cpp
   A module of spcolo
   C interface to call Cholesky factorization

   October 2010, Kenta KATO
   ************************************************************ */

#include <iostream>
#include <cstdio>
#include "spcoloextmodule.h"
#include "spcolo_SparseMatrix.h"
#include "spcolo_util.h"


#define MX_DEBUG 0

/* ============================================================
   Cholesky factorization

   This function is destructive.
   The decomposed matrix is overwritten to A.
   ============================================================ */
SparseMatrix* spcolo_cholesky(SparseMatrix* Input)
{
    int size_n = Input->getSize(SPMATRIX_ROW);
    int T_rowind[size_n];
    double T_values[size_n];
    int T_numNonzero;
    int start_ptr, end_ptr;
    double val;
    double pivotVal;
    int size_alloc;
    double ret_array[size_n];
    SparseMatrix* newMatrix = NULL;
    double tmp_val;
    int tmp_ind;
    double zeroValue = 1.0E-8;
    int rowPosition[size_n];
    for (int index = 0; index < size_n; index++) {
        rowPosition[index] = index;
    }

    SparseMatrix* A = new SparseMatrix(size_n, size_n, Input->getNumNonzero());
    int input_row, input_col;
    double input_val;
    Input->resetIterator();
    while (Input->getNext(&input_row, &input_col, &input_val)) {
        A->pushBack(input_row, input_col, input_val);
    }

    /* Start decomposition */
    for (int index = 0; index < size_n; index++) {
        T_numNonzero = 0;
        start_ptr = A->colptr[index];
        end_ptr = A->colptr[index + 1];
        pivotVal = 0;
        for (int ptr = start_ptr; ptr < end_ptr; ptr++) {
            if (A->rowind[ptr] == index) {
                pivotVal = A->values[ptr];
                break;
            }
        }

        for (int ptr = start_ptr; ptr < end_ptr; ptr++) {
            /* Make T */
            if (A->rowind[ptr] < index) {
                continue;
            }

            if (A->rowind[ptr] == index) {
                val = (1 / pivotVal) - 1;
            } else {
                val = -1 * A->values[ptr] / pivotVal;
            }

            if (abs(val) > zeroValue) {
                T_rowind[T_numNonzero] = A->rowind[ptr];
                T_values[T_numNonzero] = val;
                T_numNonzero++;
            }
        }

        /* Update matrix A */
        size_alloc = max(size_n * (size_n + 1) / 2, A->getNumNonzero());
        newMatrix = new SparseMatrix(size_n, size_n, size_alloc);

        for (int col = 0; col < size_n; col++) {
            start_ptr = A->colptr[col];
            end_ptr = A->colptr[col + 1];

            /* Initialize buffer */
            for (int row = 0; row < size_n; row++) {
                ret_array[row] = 0;
            }

            /* Calculation T * A */
            for (int ptr_target = start_ptr; ptr_target < end_ptr; ptr_target++) {
                ret_array[A->rowind[ptr_target]] += A->values[ptr_target];
                if (A->rowind[ptr_target] == index) {
                    val = A->values[ptr_target];
                    for (int ptr_T = 0; ptr_T < T_numNonzero; ptr_T++) {
                        ret_array[T_rowind[ptr_T]] += val * T_values[ptr_T];
                    }
                }
            }

            /* Copy ret_array to newMatrix */
            for (int row = 0; row < size_n; row++) {
                //if (abs(ret_array[row]) > zeroValue) {
                    newMatrix->pushBack(row, col, ret_array[row]);
                    //}
            }
        }

        delete A;
        A = newMatrix;
    }

#if MX_DEBUG
    rMessage("");
#endif

    return A;
}
