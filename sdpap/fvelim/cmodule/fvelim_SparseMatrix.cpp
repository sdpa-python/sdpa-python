/*
C class definition of `SparseMatrix`
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

TODO: Merge with `SparseMatrix` contained in SPCOLO Extension folder
*/

#include <cstdio>
#include <new>
#include "fvelim_SparseMatrix.h"


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



/* ********************************************************************************
   Matrix Operation
   ******************************************************************************** */

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define max(x, y) ((x) >= (y) ? (x) : (y))
#define min(x, y) ((x) <= (y) ? (x) : (y))

/* ============================================================
   Sub Matrix
   ============================================================ */
SparseMatrix* Matrix_sub(SparseMatrix* L, SparseMatrix* R, double zeroValue)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if ((L->getSize(SPMATRIX_ROW) != R->getSize(SPMATRIX_ROW))
        || (L->getSize(SPMATRIX_COL) != R->getSize(SPMATRIX_COL))) {
        rError("Matrix_sub(): Matrix size is not match");
    }

    int size_row = L->getSize(SPMATRIX_ROW);
    int size_col = L->getSize(SPMATRIX_COL);
    int size_alloc = max(L->getNumNonzero(), R->getNumNonzero());

    SparseMatrix* ret = new SparseMatrix(size_row, size_col, size_alloc);

    bool isEnd_L, isEnd_R;
    int col_L, col_R;
    int row_L, row_R;
    double val_L, val_R;

    L->resetIterator();
    R->resetIterator();

    isEnd_L = L->getNext(&row_L, &col_L, &val_L);
    isEnd_R = R->getNext(&row_R, &col_R, &val_R);

    int pos_L = col_L * size_row + row_L;
    int pos_R = col_R * size_row + row_R;
    double setVal;

    while (isEnd_L == true && isEnd_R == true) {
        if (pos_L == pos_R) {
            setVal = val_L - val_R;
            if (abs(setVal) > zeroValue) {
                ret->pushBack(row_L, col_L, setVal);
            }
            isEnd_L = L->getNext(&row_L, &col_L, &val_L);
            isEnd_R = R->getNext(&row_R, &col_R, &val_R);
            pos_L = col_L * size_row + row_L;
            pos_R = col_R * size_row + row_R;
        } else if (pos_L < pos_R) {
            ret->pushBack(row_L, col_L, val_L);
            isEnd_L = L->getNext(&row_L, &col_L, &val_L);
            pos_L = col_L * size_row + row_L;
        } else {
            ret->pushBack(row_R, col_R, -val_R);
            isEnd_R = R->getNext(&row_R, &col_R, &val_R);
            pos_R = col_R * size_row + row_R;
        }
    }

    while (isEnd_L == true) {
        ret->pushBack(row_L, col_L, val_L);
        isEnd_L = L->getNext(&row_L, &col_L, &val_L);
    }

    while (isEnd_R == true) {
        ret->pushBack(row_R, col_R, -val_R);
        isEnd_R = R->getNext(&row_R, &col_R, &val_R);
    }



    return ret;
}


/* ============================================================
   Multiply matrix * vector
   ============================================================ */
SparseMatrix* Matrix_mulMV(SparseMatrix* matL, SparseMatrix* vecR, double zeroValue)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (vecR->getSize(SPMATRIX_COL) != 1) {
        rError("Matrix_mulMV(): vecR must be vector");
    }
    if (matL->getSize(SPMATRIX_COL) != vecR->getSize(SPMATRIX_ROW)) {
        rError("Matrix_mulMV(): size invalid");
    }

    /* --------------------------------------------------
       make buffer to store result
       -------------------------------------------------- */
    int size_ret = matL->getSize(SPMATRIX_ROW);
    double ret_array[size_ret];
    for (int i = 0; i < size_ret; i++) {
        ret_array[i] = 0;
    }

    /* --------------------------------------------------
       multiply and store result to ret_array
       -------------------------------------------------- */
    int get_ind, set_ind;
    SparseMatrix* vecL;
    int numNonzero_R = vecR->getNumNonzero();
    int start_L, end_L;
    double valL, valR;
    for (int i = 0; i < numNonzero_R; i++) {
        get_ind = vecR->rowind[i];
        valR = vecR->values[i];
        start_L = matL->colptr[get_ind];
        end_L = matL->colptr[get_ind + 1];
        for (int j = start_L; j < end_L; j++) {
            set_ind = matL->rowind[j];
            valL = matL->values[j];
            ret_array[set_ind] += valL * valR;
        }
    }

    /* --------------------------------------------------
       copy ret_array to SparseMatrix
       -------------------------------------------------- */
    int numNonzero_ret = 0;
    for (int i = 0; i < size_ret; i++) {
        if (ret_array[i] != 0) {
            numNonzero_ret++;
        }
    }

    SparseMatrix* ret = new SparseMatrix(size_ret, 1, numNonzero_ret);
    for (int i = 0; i < size_ret; i++) {
        if (abs(ret_array[i]) > zeroValue) {
            ret->pushBack(i, 0, ret_array[i]);
        }
    }

    return ret;
}


/* ============================================================
   Multiply vector * matrix
   ============================================================ */
SparseMatrix* Matrix_mulVM(SparseMatrix* vecL, SparseMatrix* matR, double zeroValue)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (vecL->getSize(SPMATRIX_COL) != 1) {
        rError("Matrix_mulVM(): vecL must be vector");
    }
    if (matR->getSize(SPMATRIX_ROW) != vecL->getSize(SPMATRIX_ROW)) {
        rError("Matrix_mulVM(): size invalid");
    }
    /* -------------------------------------------------- */

    int size_ret = matR->getSize(SPMATRIX_COL);
    int numNonzero_L = vecL->getNumNonzero();

    SparseMatrix* ret = new SparseMatrix(size_ret, 1, numNonzero_L);

    int start_R, end_R;
    int ptr_L, ptr_R;
    double product;
    for (int i = 0; i < size_ret; i++) {
        start_R = matR->colptr[i];
        end_R = matR->colptr[i + 1];
        ptr_L = 0;
        ptr_R = start_R;
        product = 0;
        while (ptr_L < numNonzero_L && ptr_R < end_R) {
            if (vecL->rowind[ptr_L] == matR->rowind[ptr_R]) {
                product += vecL->values[ptr_L] * matR->values[ptr_R];
                ptr_L++;
                ptr_R++;
            } else if (vecL->rowind[ptr_L] < matR->rowind[ptr_R]) {
                ptr_L++;
            } else {
                ptr_R++;
            }
        }
        if (abs(product) > zeroValue) {
            ret->pushBack(i, 0, product);
        }
    }

    return ret;
}


/* ============================================================
   Multiply Matrix * Matrix
   ============================================================ */
SparseMatrix* Matrix_mulMM(SparseMatrix* matL, SparseMatrix* matR, double zeroValue)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (matL->getSize(SPMATRIX_COL) != matR->getSize(SPMATRIX_ROW)) {
        rError("Matrix_mulMM(): size invalid");
    }

    /* --------------------------------------------------
       make buffer to store result
       -------------------------------------------------- */
    int size_row = matL->getSize(SPMATRIX_ROW);
    double ret_array[size_row];

    int size_col = matR->getSize(SPMATRIX_COL);
    int size_alloc = matL->getNumNonzero() + matR->getNumNonzero();

    SparseMatrix* ret = new SparseMatrix(size_row, size_col, size_alloc);

    /* --------------------------------------------------
       make col_vector to ret_array and copy it to SparseMatrix
       -------------------------------------------------- */
    int start_R, end_R;
    int start_L, end_L;
    int get_ind, set_ind;
    double valR, valL;
    for (int col = 0; col < size_col; col++) {
        start_R = matR->colptr[col];
        end_R = matR->colptr[col + 1];

        /* ----------------------------------------
           initialize buffer
           ---------------------------------------- */
        for (int row = 0; row < size_row; row++) {
            ret_array[row] = 0;
        }

        /* ----------------------------------------
           multiply and store result to ret_array
           ---------------------------------------- */
        for (int i = start_R; i < end_R; i++) {
            get_ind = matR->rowind[i];
            valR = matR->values[i];
            start_L = matL->colptr[get_ind];
            end_L = matL->colptr[get_ind + 1];
            for (int j = start_L; j < end_L; j++) {
                set_ind = matL->rowind[j];
                valL = matL->values[j];
                ret_array[set_ind] += valL * valR;
            }
        }

        /* ----------------------------------------
           copy ret_array to SparseMatrix
           ---------------------------------------- */
        for (int row = 0; row < size_row; row++) {
            if (abs(ret_array[row]) > zeroValue) {
                ret->pushBack(row, col, ret_array[row]);
            }
        }

    }

    return ret;
}




/* ============================================================
   Inner Product of two vectors
   ============================================================ */
double Matrix_innerProduct(SparseMatrix* vec1, SparseMatrix* vec2, double zeroValue)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (vec1->getSize(SPMATRIX_COL) != 1 || vec2->getSize(SPMATRIX_COL) != 1) {
        rError("Matrix_innerProduct(): vec1 & vec2 must be vector");
    }
    if (vec1->getSize(SPMATRIX_ROW) != vec2->getSize(SPMATRIX_ROW)) {
        rError("Matrix_innerProduct(); size invalid");
    }
    /* -------------------------------------------------- */

    int numNonzero_1 = vec1->getNumNonzero();
    int numNonzero_2 = vec1->getNumNonzero();
    int ptr_1 = 0;
    int ptr_2 = 0;
    double ret = 0;
    while (ptr_1 < numNonzero_1 || ptr_2 < numNonzero_2) {
        if (vec1->rowind[ptr_1] == vec2->rowind[ptr_2]) {
            ret += vec1->values[ptr_1] * vec2->values[ptr_2];
            ptr_1++;
            ptr_2++;
        } else if (vec1->rowind[ptr_1] < vec2->rowind[ptr_2]) {
            ptr_1++;
        } else {
            ptr_2++;
        }
    }

    if (abs(ret) > zeroValue) {
        return ret;
    } else {
        return 0;
    }
}


/* ============================================================
   Calc (UnitUpperTrianglarMatrix)^(-1) * (Vector) by solve equation
   ============================================================ */
SparseMatrix* Matrix_solveEquationMV(SparseMatrix* mat, SparseMatrix* vec, double zeroValue)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (vec->getSize(SPMATRIX_COL) != 1) {
        rError("Matrix_solveEquationMV(): vec must be vector");
    }

    if (mat->getSize(SPMATRIX_ROW) != mat->getSize(SPMATRIX_COL)){
        rError("Matrix_solveEquationMV(): mat must be square");
    }

    if (mat->getSize(SPMATRIX_COL) != vec->getSize(SPMATRIX_ROW)){
        rError("Matrix_solveEquationMV(): size invalid");
    }


    /* --------------------------------------------------
       make buffer and copy vec to rem
       -------------------------------------------------- */
    int size_ret = vec->getSize(SPMATRIX_ROW);
    int numNonzero_vec = vec->getNumNonzero();
    double rem[size_ret];
    for (int row = 0; row < size_ret; row++) {
        rem[row] = 0;
    }

    for (int i = 0; i < numNonzero_vec; i++) {
        rem[vec->rowind[i]] = vec->values[i];
    }

    /* --------------------------------------------------
       solve equation
       -------------------------------------------------- */
    int start_mat, end_mat;
    for (int index = size_ret - 1; index >= 0; index--) {
        start_mat = mat->colptr[index];
        end_mat = mat->colptr[index + 1];
        for (int ptr = start_mat; ptr < end_mat; ptr++) {
            if (mat->rowind[ptr] == index) {
                break;
            }
            rem[mat->rowind[ptr]] -= rem[index] * mat->values[ptr];
        }
    }

    /* --------------------------------------------------
       copy rem to SparseMatrix
       -------------------------------------------------- */
    int numNonzero_ret = 0;
    for (int row = 0; row < size_ret; row++) {
        if (abs(rem[row]) > zeroValue) {
            numNonzero_ret++;
        }
    }

    SparseMatrix* ret = new SparseMatrix(size_ret, 1, numNonzero_ret);

    for (int row = 0; row < size_ret; row++) {
        if (abs(rem[row]) > zeroValue) {
            ret->pushBack(row, 0, rem[row]);
        }
    }

    return ret;
}


/* ============================================================
   Calc (Vector)^T * (UpperTrianglarMatrix)^(-1) by solve equation
   ============================================================ */
SparseMatrix* Matrix_solveEquationVM(SparseMatrix* vec, SparseMatrix* mat, double zeroValue)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (vec->getSize(SPMATRIX_COL) != 1) {
        rError("Matrix_solveEquationMV(): vec must be vector");
    }

    if (mat->getSize(SPMATRIX_ROW) != mat->getSize(SPMATRIX_COL)){
        rError("Matrix_solveEquationMV(): mat must be square");
    }

    if (mat->getSize(SPMATRIX_COL) != vec->getSize(SPMATRIX_ROW)){
        rError("Matrix_solveEquationMV(): size invalid");
    }

    /* --------------------------------------------------
       make buffer to store result
       -------------------------------------------------- */
    int size_ret = vec->getSize(SPMATRIX_ROW);
    int numNonzero_vec = vec->getNumNonzero();
    double rem[size_ret];
    for (int row = 0; row < size_ret; row++) {
        rem[row] = 0;
    }

    for (int i = 0; i < numNonzero_vec; i++) {
        rem[vec->rowind[i]] = vec->values[i];
    }

    /* --------------------------------------------------
       solve equation
       -------------------------------------------------- */
    int start_mat, end_mat;
    for (int index = 0; index < size_ret; index++) {
        start_mat = mat->colptr[index];
        end_mat = mat->colptr[index + 1];
        for (int ptr = start_mat; ptr < end_mat; ptr++) {
            if (mat->rowind[ptr] == index) {
                rem[index] /= mat->values[ptr];
                break;
            }
            rem[index] -= rem[mat->rowind[ptr]] * mat->values[ptr];
        }
    }

    /* --------------------------------------------------
       copy rem to SparseMatrix
       -------------------------------------------------- */
    int numNonzero_ret = 0;
    for (int row = 0; row < size_ret; row++) {
        if (rem[row] != 0) {
            numNonzero_ret++;
        }
    }

    SparseMatrix* ret = new SparseMatrix(size_ret, 1, numNonzero_ret);

    for (int row = 0; row < size_ret; row++) {
        if (abs(rem[row]) > zeroValue) {
            ret->pushBack(row, 0, rem[row]);
        }
    }

    return ret;
}

