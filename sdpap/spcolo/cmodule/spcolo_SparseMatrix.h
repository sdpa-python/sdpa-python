/* ********************************************************************************
   spcolo_SparseMatrix.h
   A module of spcolo
   C interface to call free variable elimination

   December 2010, Kenta KATO
   ******************************************************************************** */

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
