/* ********************************************************************************
   spcolo_SparseMatrix.cpp
   A module of spcolo
   C interface to call Cholesky factorization

   December 2010, Kenta KATO
   ******************************************************************************** */

#include <cstdio>
#include <new>
#include "spcolo_SparseMatrix.h"


/* ************************************************************
   SparseMatrix
   ************************************************************ */

/* ==================================================
   Constructor
   ================================================== */
SparseMatrix::SparseMatrix(int size_row, int size_col, int size_alloc) :
    size_row(size_row), size_col(size_col), size_alloc(size_alloc)
{
    NewArray(colptr, int, size_col + 1);
    /*
    colptr = new int[size_col + 1];
    */
    for(int i = 0; i <= size_col; i++){
        colptr[i] = 0;
    }

    int allocation = size_alloc;
    if (size_alloc == 0) {
        allocation = 1;
    }
    NewArray(values, double, allocation);
    NewArray(rowind, int, allocation);

    return;
}


/* ==================================================
   getValue
   ================================================== */
double SparseMatrix::getValue(int row, int col)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (row >= size_row || col >= size_col) {
        rError("SparseMatrix::getValue(): size over");
    }

    /* ------------------------------
       Binary search
       ------------------------------ */
    int start;
    int end;
    int pos;

    start = colptr[col];
    end = colptr[col + 1];

    if (start == end) {
        return 0;
    }

    pos = (start + end) / 2;

    while (start != end) {
        if (rowind[pos] == row) {
            return values[pos];
        }

        if (row < rowind[start] || rowind[end - 1] < row) {
            return 0;
        }

        if (rowind[pos] < row) {
            start = pos + 1;
        }else{
            end = pos - 1;
        }
        pos = (start + end) / 2;
    }

    if (rowind[pos] == row) {
        return values[pos];
    }

    return 0;
}


/* ============================================================
   getVector
   ============================================================ */
SparseMatrix* SparseMatrix::getVector(int col)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (col >= size_col) {
        rError("SparseMatrix::getVector(): size over");
    }

    /* --------------------------------------------------
       get vector
       -------------------------------------------------- */
    int size_alloc = getNumNonzero(col);
    if (size_alloc == 0) {
        size_alloc = 1;
    }
    SparseMatrix* ret = new SparseMatrix(size_row, 1, size_alloc);

    int start = colptr[col];
    int end = colptr[col + 1];

    if (start == end) {
        return ret;
    }

    for (int i = start; i < end; i++) {
        ret->values[i - start] = this->values[i];
        ret->rowind[i - start] = this->rowind[i];
    }

    ret->colptr[0] = 0;
    ret->colptr[1] = end - start;

    return ret;
}


/* ============================================================
   Get next value by Iterator
   pos_* is indexed 1,2,...(not 0,1,...)
   [return]
   false: not have element
   ============================================================ */
bool SparseMatrix::getNext(int* row, int* col, double* val)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (iter == getNumNonzero()) {
        *row = size_row;
        *col = size_col;
        *val = 0;
        return false;
    }

    *row = rowind[iter];
    *col = cur_ptr;
    *val = values[iter];

    iter++;
    if (cur_ptr < size_col) {
        while (iter >= colptr[cur_ptr + 1]) {
            cur_ptr++;
            if (cur_ptr == size_col) {
                break;
            }
            if (cur_ptr + 1 >= size_col + 1) {
                printf("size_col = %d\n", size_col);
                printf("memory leak!: %d\n", cur_ptr);
            }
        }
    }

    return true;
}


/* ============================================================
   pushElement
   ============================================================ */
void SparseMatrix::pushBack(int row, int col, double val)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (row >= size_row || col >= size_col) {
        rError("SparseMatrix::pushBack(): size over");
    }

    if (val == 0) {
        return;
    }

    int pushPoint = getNumNonzero();

    if (pushPoint == size_alloc) {
        /* --------------------------------------------------
           Reallocation memory
           -------------------------------------------------- */
        double* new_values;
        int* new_rowind;
        NewArray(new_values, double, size_alloc * 2);
        NewArray(new_rowind, int, size_alloc * 2);
        for (int i = 0; i < pushPoint; i++) {
            new_values[i] = values[i];
            new_rowind[i] = rowind[i];
        }
        DeleteArray(values);
        DeleteArray(rowind);
        values = new_values;
        rowind = new_rowind;
        size_alloc *= 2;
    }

    /* --------------------------------------------------
       Push new element
       -------------------------------------------------- */
    values[pushPoint] = val;
    rowind[pushPoint] = row;
    for (int i = col + 1; i <= size_col; i++) {
        colptr[i]++;
    }

    return;
}


/* ============================================================
   setVector
   ============================================================ */
void SparseMatrix::setVector(int col, SparseMatrix* vec)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (vec->getSize(SPMATRIX_COL) != 1) {
        rError("SparseMatrix::setVector(): vector size_col must be 1");
    }
    if (vec->getSize(SPMATRIX_ROW) != size_row) {
        rError("SparseMatrix::setVector(): vector size is invalid.");
    }

    /* --------------------------------------------------
       Set vector
       -------------------------------------------------- */
    int* old_colptr = colptr;
    int* old_rowind = rowind;
    double* old_values = values;

    NewArray(colptr, int, size_col + 1);
    NewArray(rowind, int, size_alloc);
    NewArray(values, double, size_alloc);

    for (int i = 0; i <= size_col; i++) {
        colptr[i] = 0;
    }

    int set_ptr = 0;

    /* --------------------------------------------------
       Copy former matrix
       -------------------------------------------------- */
    for (int i = 0; i < old_colptr[col]; i++) {
        while (i >= old_colptr[set_ptr + 1]) {
            set_ptr++;
        }
        pushBack(old_rowind[i], set_ptr, old_values[i]);
    }

    /* --------------------------------------------------
       Insert new vector
       -------------------------------------------------- */
    int numNonzero_vec = vec->getNumNonzero();
    for (int i = 0; i < numNonzero_vec; i++) {
        pushBack(vec->rowind[i], col, vec->values[i]);
    }

    /* --------------------------------------------------
       Copy latter matrix
       -------------------------------------------------- */
    int numNonzero = old_colptr[size_col];
    for (int i = old_colptr[col + 1]; i < numNonzero; i++) {
        while (i >= old_colptr[set_ptr + 1]) {
            set_ptr++;
        }
        pushBack(old_rowind[i], set_ptr, old_values[i]);
    }

    DeleteArray(old_values);
    DeleteArray(old_rowind);
    DeleteArray(old_colptr);

    return;
}



/* ============================================================
   permute 2 elements in vector
   ============================================================ */
void SparseMatrix::permuteVector(int index1, int index2)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (size_col != 1) {
        rError("SparseMatrix::permuteVector(): the Matrix type must be vector");
    }

    if (index1 >= size_row || index2 >= size_row) {
        rError("SparseMatrix::permuteVector(): size over");
    }

    /* --------------------------------------------------
       permute
       -------------------------------------------------- */
    if (index1 == index2) {
        return;
    }

    int indexF, indexL;
    if (index1 < index2) {
        indexF = index1;
        indexL = index2;
    } else {
        indexF = index2;
        indexL = index1;
    }

    double valL = getValue(indexL, 0);

    int* old_rowind = rowind;
    double* old_values = values;

    int numNonzero = getNumNonzero();

    NewArray(rowind, int, size_alloc);
    NewArray(values, double, size_alloc);

    int index_old = 0;
    int index_new = 0;
    while (old_rowind[index_old] < indexF && index_old < numNonzero) {
        rowind[index_new] = old_rowind[index_old];
        values[index_new] = old_values[index_old];
        index_old++;
        index_new++;
    }

    if (index_old == numNonzero) {
        /* Both of permutation value are 0 */
        return;
    }

    double valF;
    if (old_rowind[index_old] == indexF) {
        valF = old_values[index_old];
        if (valL != 0) {
            rowind[index_new] = indexF;
            values[index_new] = valL;
            index_new++;
        }
        index_old++;
    } else {
        valF = 0;
        if (valL != 0) {
            rowind[index_new] = indexF;
            values[index_new] = valL;
            index_new++;
        }
    }

    while (old_rowind[index_old] < indexL && index_old < numNonzero) {
        rowind[index_new] = old_rowind[index_old];
        values[index_new] = old_values[index_old];
        index_old++;
        index_new++;
    }

    if (index_old < numNonzero && old_rowind[index_old] == indexL) {
        index_old++;
    }

    if (valF != 0) {
        rowind[index_new] = indexL;
        values[index_new] = valF;
        index_new++;
    }

    while (index_old < numNonzero) {
        rowind[index_new] = old_rowind[index_old];
        values[index_new] = old_values[index_old];
        index_old++;
        index_new++;
    }

    DeleteArray(old_rowind);
    DeleteArray(old_values);

    return;
}



/* ============================================================
   permute 2 vectors in matrix
   ============================================================ */
void SparseMatrix::permuteMatrix(int index1, int index2)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (index1 >= size_col || index2 >= size_col) {
        rError( "SparseMatrix::permute(): size over" );
    }

    /* --------------------------------------------------
       permute
       -------------------------------------------------- */
    SparseMatrix* vec1 = getVector(index1);
    SparseMatrix* vec2 = getVector(index2);
    setVector(index2, vec1);
    setVector(index1, vec2);
    delete vec1;
    delete vec2;

    return;
}


/* ============================================================
   print
   ============================================================ */
void SparseMatrix::print()
{
    int index_ind = 0;
    int row, col;
    for (int col = 0; col < size_col; col++) {
        int numNonzero = getNumNonzero(col);
        for (int count = 0; count < numNonzero; count++) {
            row = rowind[index_ind];
            printf("(%d,%d) = %f\n", row, col, values[index_ind]);
            index_ind++;
        }
    }

    printf("--------------------\n");

    return;
}


/* ============================================================
   Destructor
   ============================================================ */
SparseMatrix::~SparseMatrix()
{
    DeleteArray(values);
    DeleteArray(rowind);
    DeleteArray(colptr);
    return;
}
