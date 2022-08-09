/*
C class definition of `SparseMatrix` and some friend functions
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
*/

#ifndef __SPCOLO_SPARSEMATRIX_H__
#define __SPCOLO_SPARSEMATRIX_H__

#include "spcolo_util.h"
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

/* ============================================================
   SparseMatrix
   Sparse Matrix with Compressed Column Storage format (CCS)
   values: value of nonzero elements
   rowind: index of row of nonzero elements
   colptr: the first point of each column
   index is zero base (0, 1, 2, ...)
   ============================================================ */

typedef enum {
    SPMATRIX_ROW,
    SPMATRIX_COL
}VectorType;


class SparseMatrix
{
 public:
    /* Matrix Size */
    int size_row;
    int size_col;

    /* Matrix Elements */
    double* values;
    int* rowind;
    int* colptr;

    /* Allocated memory infomation */
    int size_alloc;

    /* Iterator */
    int iter;
    int cur_ptr;

    /* --------------------------------------------------
       Constructor
       -------------------------------------------------- */
    SparseMatrix(int size_row, int size_col, int size_alloc);

    /* --------------------------------------------------
       Reset Iterator
       -------------------------------------------------- */
    void resetIterator()
    {
        iter = 0;
        cur_ptr = 0;
    }

    /* --------------------------------------------------
       Getter functions
       -------------------------------------------------- */
    double getValue(int row, int col);
    SparseMatrix* getVector(int index);

    int getSize(VectorType type)
    {
        return (type == SPMATRIX_ROW) ? size_row : size_col;
    }

    /* Get number of nonzero elements in vector */
    int getNumNonzero(int index)
    {
        return colptr[index + 1] - colptr[index];
    }

    /* Get number of nonzero elements in matrix */
    int getNumNonzero()
    {
        return colptr[size_col];
    }

    /* --------------------------------------------------
       Get next value by Iterator
       [return]
       false: not have element
       -------------------------------------------------- */
    bool getNext(int* row, int* col, double* val);

    /* --------------------------------------------------
       Setter functions
       -------------------------------------------------- */
    void setVector(int index, SparseMatrix* vector);

    void pushBack(int row, int col, double val);

    /* --------------------------------------------------
       permute 2 elements in vector
       -------------------------------------------------- */
    void permuteVector(int index1, int index2);

    /* --------------------------------------------------
       permute 2 vectors in matrix
       -------------------------------------------------- */
    void permuteMatrix(int index1, int index2);

    /* --------------------------------------------------
       print
       Print nonzero elements
       -------------------------------------------------- */
    void print();

    /* --------------------------------------------------
       Destructor
       -------------------------------------------------- */
    ~SparseMatrix();
};

#endif
