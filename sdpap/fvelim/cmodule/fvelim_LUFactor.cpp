/*
C class definition of `LUFactor` (class to perform LU factorization of a `SparseMatrix`)
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

#include "fvelim_LUFactor.h"
#include <cstdio>
#include <limits.h>

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define max(x, y) ((x) >= (y) ? (x) : (y))

#define LU_DEBUG 0
#define LU_TRACE 0

/* ============================================================
   Constructor
   ============================================================ */
LUFactor::LUFactor(SparseMatrix* A_f, SparseMatrix* A_lqs, double rho, double zeroValue):
    size_Kf(A_f->getSize(SPMATRIX_COL)), size_Alqs(A_lqs->getSize(SPMATRIX_COL)), size_m(A_f->getSize(SPMATRIX_ROW)), rho(rho), zeroValue(zeroValue)
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (rho < 0 || rho > 1) {
        rError("LUInfo::LUInfo(): rho must be 0 < rho < 1");
    }
    /* -------------------------------------------------- */

#if LU_TRACE
    rMessage("");
#endif

    int row, col;
    double val;

    /* --------------------------------------------------
       Initialize pivot_Af
       -------------------------------------------------- */
    pivot_Af = new SparseMatrix(size_m, size_Kf, A_f->getNumNonzero());
    A_f->resetIterator();
    while (A_f->getNext(&row, &col, &val) != false) {
        pivot_Af->pushBack(row, col, val);
    }

#if LU_TRACE
    rMessage("");
#endif

    /* --------------------------------------------------
       Initialize pos_Annz and num_Annz
       -------------------------------------------------- */
    NewArray(pos_Annz, int*, size_m);
    NewArray(num_Annz, int, size_m);
    for (int i = 0; i < size_m; i++) {
        NewArray(pos_Annz[i], int, size_Alqs + 1);
        num_Annz[i] = 0;
    }

    A_lqs->resetIterator();
    while (A_lqs->getNext(&row, &col, &val) != false) {
        pos_Annz[row][num_Annz[row]] = col;
        num_Annz[row]++;
    }

    for (int i = 0; i < size_m; i++) {
        pos_Annz[i][num_Annz[i]] = -1;
    }

#if LU_TRACE
    rMessage("");
    printf("size_m = %d\n", size_m);
#endif


    /* --------------------------------------------------
       Initialize LInvP
       -------------------------------------------------- */
    LInvP = new SparseMatrix(size_m, size_m, 2 * size_m);

#if LU_TRACE
    rMessage("");
#endif

    /* Make Identity Matrix */
    for (int i = 0; i < size_m; i++) {
        LInvP->pushBack(i, i, 1);
    }

#if LU_TRACE
    rMessage("");
#endif

    /* --------------------------------------------------
       Initialize row, col position
       -------------------------------------------------- */
    NewArray(colPosition, int, size_Kf);

    for (int i = 0; i < size_Kf; i++) {
        colPosition[i] = i;
    }

#if LU_TRACE
    rMessage("");
#endif

    /* --------------------------------------------------
       Initialize P and Q
       -------------------------------------------------- */
    NewArray(P, int, size_Kf);
    NewArray(Q, int, size_Kf);
    P[0] = -1;
    Q[0] = -1;

    rank_Af = 0;

    return;
}


/* ============================================================
   setPivot
   ============================================================ */
void LUFactor::setPivot(int pivotRow, int pivotCol, double pivotVal, int index)
{
    int T_rowind[pivot_Af->getNumNonzero(pivotCol)];
    double T_values[pivot_Af->getNumNonzero(pivotCol)];
    int T_numNonzero = 0;
    int start_ptr;
    int end_ptr;
    double val;

    /* --------------------------------------------------
       make T
       -------------------------------------------------- */
    start_ptr = pivot_Af->colptr[pivotCol];
    end_ptr = pivot_Af->colptr[pivotCol + 1];

    for (int ptr = start_ptr; ptr < end_ptr; ptr++) {
        if (pivot_Af->rowind[ptr] < index) {
            continue;
        }

        if (pivot_Af->rowind[ptr] == pivotRow) {
            val = (1 / pivotVal) - 1;
        } else {
            val = -1 * pivot_Af->values[ptr] / pivotVal;
        }

        if (abs(val) > zeroValue) {
            T_rowind[T_numNonzero] = pivot_Af->rowind[ptr];
            T_values[T_numNonzero] = val;
            T_numNonzero++;
        }
    }


    /* --------------------------------------------------
       set P, Q
       -------------------------------------------------- */
    P[index] = pivotRow;
    if (index != size_Kf - 1) {
        P[index + 1] = -1;
    }

    bool isSetQ = false;
    for (int i = 0; i < size_Kf; i++) {
        if (colPosition[i] == pivotCol) {
            Q[index] = i;
            if (index != size_Kf - 1) {
                Q[index + 1] = -1;
            }
            int tmp = colPosition[index];
            colPosition[index] = colPosition[i];
            colPosition[i] = tmp;
            isSetQ = true;
            break;
        }
    }
    if (isSetQ == false) {
        rError("LUInfo::setPivot(): an error occured");
    }

    /* --------------------------------------------------
       update pos_Annz
       -------------------------------------------------- */
    int* new_pos;
    int* pos1;
    int* pos2;
    int ptr1, ptr2;
    int ptr_new;

    for (int ptr = start_ptr; ptr < end_ptr; ptr++) {
        if (pivot_Af->rowind[ptr] == pivotRow) {
            continue;
        }

        NewArray(new_pos, int, size_Alqs + 1);
        pos1 = pos_Annz[pivotRow];
        pos2 = pos_Annz[pivot_Af->rowind[ptr]];
        ptr1 = 0;
        ptr2 = 0;
        ptr_new = 0;
        while (pos1[ptr1] != -1 || pos2[ptr2] != -1) {
            if (pos1[ptr1] != -1 && pos2[ptr2] != -1 && pos1[ptr1] == pos2[ptr2]) {
                new_pos[ptr_new] = pos1[ptr1];
                ptr1++;
                ptr2++;
                ptr_new++;
            } else if (pos1[ptr1] != -1 && (pos2[ptr2] == -1 || pos1[ptr1] < pos2[ptr2])) {
                new_pos[ptr_new] = pos1[ptr1];
                ptr1++;
                ptr_new++;
            } else {
                new_pos[ptr_new] = pos2[ptr2];
                ptr2++;
                ptr_new++;
            }
        }
        new_pos[ptr_new] = -1;
        DeleteArray(pos_Annz[pivot_Af->rowind[ptr]]);
        pos_Annz[pivot_Af->rowind[ptr]] = new_pos;
        num_Annz[pivot_Af->rowind[ptr]] = ptr_new;
    }

    {
        int* tmp_pos = pos_Annz[index];
        int tmp_nnz = num_Annz[index];
        pos_Annz[index] = pos_Annz[pivotRow];
        num_Annz[index] = num_Annz[pivotRow];
        pos_Annz[pivotRow] = tmp_pos;
        num_Annz[pivotRow] = tmp_nnz;
    }

    /* --------------------------------------------------
       Update pivot_Af and LInvP by multiplying T and P
       -------------------------------------------------- */
    /* buffer to store result */
    double ret_array[size_m];

    int size_alloc;
    SparseMatrix* newMatrix;

    double tmp_val;

    /* ----------------------------------------
       update pivot_Af
       ---------------------------------------- */
    size_alloc = 2 * pivot_Af->getNumNonzero();

    newMatrix = new SparseMatrix(size_m, size_Kf, size_alloc);

    for (int col = 0; col < size_Kf; col++) {
        start_ptr = pivot_Af->colptr[col];
        end_ptr = pivot_Af->colptr[col + 1];

        /*
         * initialize buffer
         */
        for (int row = 0; row < size_m; row++) {
            ret_array[row] = 0;
        }

        /*
         * calculation T * pivot_Af
         */
        for (int ptr_target = start_ptr; ptr_target < end_ptr; ptr_target++) {
            ret_array[pivot_Af->rowind[ptr_target]] += pivot_Af->values[ptr_target];
            if (pivot_Af->rowind[ptr_target] == pivotRow) {
                val = pivot_Af->values[ptr_target];
                for (int ptr_T = 0; ptr_T < T_numNonzero; ptr_T++) {
                    ret_array[T_rowind[ptr_T]] += val * T_values[ptr_T];
                }
            }
        }

        /*
         * permutation
         */
        tmp_val = ret_array[index];
        ret_array[index] = ret_array[pivotRow];
        ret_array[pivotRow] = tmp_val;

        /*
         *  copy ret_array to newMatrix
         */
        for (int row = 0; row < size_m; row++) {
            if (abs(ret_array[row]) > zeroValue) {
                newMatrix->pushBack(row, col, ret_array[row]);
            }
        }
    }

    delete pivot_Af;
    pivot_Af = NULL;
    pivot_Af = newMatrix;


    /* ----------------------------------------
       update LInvP
       ---------------------------------------- */
    size_alloc = 2 * LInvP->getNumNonzero();

    newMatrix = new SparseMatrix(size_m, size_m, size_alloc);

    for (int col = 0; col < size_m; col++) {
        start_ptr = LInvP->colptr[col];
        end_ptr = LInvP->colptr[col + 1];

        /*
         * initialize buffer
         */
        for (int row = 0; row < size_m; row++) {
            ret_array[row] = 0;
        }

        /*
         * calculation T * pivot_Af
         */
        for (int ptr_target = start_ptr; ptr_target < end_ptr; ptr_target++) {
            ret_array[LInvP->rowind[ptr_target]] += LInvP->values[ptr_target];
            if (LInvP->rowind[ptr_target] == pivotRow) {
                val = LInvP->values[ptr_target];
                for (int ptr_T = 0; ptr_T < T_numNonzero; ptr_T++) {
                    ret_array[T_rowind[ptr_T]] += val * T_values[ptr_T];
                }
            }
        }

        /*
         * permutation
         */
        tmp_val = ret_array[index];
        ret_array[index] = ret_array[pivotRow];
        ret_array[pivotRow] = tmp_val;

        /*
         *  copy ret_array to newMatrix
         */
        for (int row = 0; row < size_m; row++) {
            if (abs(ret_array[row]) > zeroValue) {
                newMatrix->pushBack(row, col, ret_array[row]);
            }
        }
    }

    delete LInvP;
    LInvP = NULL;
    LInvP = newMatrix;

    return;
}


/* ============================================================
   getU
   ============================================================ */
SparseMatrix* LUFactor::getU()
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (rank_Af == 0) {
        rError("LUInfo::getU(): LU decomposition has not finished yet.");
    }
    /* -------------------------------------------------- */

    int numNonzero = pivot_Af->getNumNonzero();
    SparseMatrix* U = new SparseMatrix(rank_Af, rank_Af, numNonzero);

    /* --------------------------------------------------
       Copy pivot_Af to U
       -------------------------------------------------- */
    int start_ptr, end_ptr;
    for (int col = 0; col < rank_Af; col++) {
        start_ptr = pivot_Af->colptr[colPosition[col]];
        end_ptr = pivot_Af->colptr[colPosition[col] + 1];
        for (int ptr = start_ptr; ptr < end_ptr; ptr++) {
            if (pivot_Af->rowind[ptr] >= rank_Af) {
                continue;
            }
            U->pushBack(pivot_Af->rowind[ptr], col, pivot_Af->values[ptr]);
        }
    }

    return U;
}


/* ============================================================
   getV
   not full rank case
   ============================================================ */
SparseMatrix* LUFactor::getV()
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (rank_Af == 0){
        rError("LUInfo::getV(): LU decomposition has not finished yet.");
    }

    if (rank_Af == size_Kf) {
        return NULL;
    }
    /* -------------------------------------------------- */

    int numNonzero = 0;
    for (int col = rank_Af; col < size_Kf; col++) {
        numNonzero += pivot_Af->getNumNonzero(colPosition[col]);
    }

    SparseMatrix* V = new SparseMatrix(rank_Af, size_Kf - rank_Af, numNonzero);

    /* --------------------------------------------------
       Copy pivot_Af to V
       -------------------------------------------------- */
    int start_ptr, end_ptr;
    for (int col = rank_Af; col < size_Kf; col++) {
        start_ptr = pivot_Af->colptr[colPosition[col]];
        end_ptr = pivot_Af->colptr[colPosition[col] + 1];
        for (int ptr = start_ptr; ptr < end_ptr; ptr++) {
            V->pushBack(pivot_Af->rowind[ptr], col - rank_Af, pivot_Af->values[ptr]);
        }
    }

    return V;
}


/* ============================================================
   Get L^(-1) * P
   ============================================================ */
SparseMatrix* LUFactor::getLInvP()
{
    SparseMatrix* ret_LInvP = new SparseMatrix(size_m, size_m, LInvP->getNumNonzero());
    int row, col;
    double val;

    /* --------------------------------------------------
       Copy LInvP
       -------------------------------------------------- */
    LInvP->resetIterator();
    while (LInvP->getNext(&row, &col, &val) != false) {
        ret_LInvP->pushBack(row, col, val);
    }

    return ret_LInvP;
}


/* ============================================================
   decompose
   LU Decomposition
   ============================================================ */
void LUFactor::decompose()
{
    /* --------------------------------------------------
       Error check
       -------------------------------------------------- */
    if (P[0] != -1) {
        rMessage("LUInfo::decompose(): Already decomposed");
        return;
    }
    /* -------------------------------------------------- */

    /* Flags whether the index was pivoted */
    bool isPivoted[size_Kf];
    for (int i = 0; i < size_Kf; i++) {
        isPivoted[i] = false;
    }

    /* Flag of independence */
    bool isIndependent = false;

#if LU_TRACE
    rMessage("");
#endif

    /* --------------------------------------------------
       Start decomposition
       -------------------------------------------------- */
    int pivotRow, pivotCol;
    double pivotVal;
    int min_num_fillin;
    int num_fillin;
    int start_ptr, end_ptr;
    int num_nnzUpper, num_nnzLower;
    double max_val;
    for (int cnt = 0; cnt < size_Kf; cnt++) {
        isIndependent = false;

        min_num_fillin = INT_MAX;
        pivotRow = -1;
        pivotCol = -1;
        pivotVal = 0;

        for (int col = 0; col < size_Kf; col++) {
            if (isPivoted[col] == true) {
                continue;
            }

            start_ptr = pivot_Af->colptr[col];
            end_ptr = pivot_Af->colptr[col + 1];

            /*
             * count number of nonzero element in col
             */
            num_nnzUpper = 0;
            for (int ptr = start_ptr; ptr < end_ptr; ptr++) {
                if (pivot_Af->rowind[ptr] >= cnt) {
                    break;
                }
                num_nnzUpper++;
            }

            num_nnzLower = pivot_Af->getNumNonzero(col) - num_nnzUpper;

            /*
             * search maximum value in col
             */
            max_val = 0;
            for (int ptr = start_ptr; ptr < end_ptr; ptr++) {
                if (pivot_Af->rowind[ptr] < cnt) {
                    continue;
                }
                if (abs(pivot_Af->values[ptr]) > max_val) {
                    max_val = abs(pivot_Af->values[ptr]);
                }
            }

            /*
             * search pivot element
             */
            for (int ptr = start_ptr; ptr < end_ptr; ptr++) {
                if (pivot_Af->rowind[ptr] < cnt) {
                    continue;
                }

                if (pivotRow == -1) {
                    // Initialize pivot
                    pivotRow = pivot_Af->rowind[start_ptr];
                    pivotCol = col;
                    pivotVal = pivot_Af->values[start_ptr];
                    isIndependent = true;
                }

                if (abs(pivot_Af->values[ptr]) >= max_val * rho) {
                    num_fillin = (num_nnzLower - 1) * num_Annz[pivot_Af->rowind[ptr]];

                    if ((num_fillin < min_num_fillin)
                        || (num_fillin == min_num_fillin
                            && abs(pivot_Af->values[ptr]) > abs(pivotVal))) {
                        pivotRow = pivot_Af->rowind[ptr];
                        pivotCol = col;
                        min_num_fillin = num_fillin;
                        pivotVal = pivot_Af->values[ptr];
                        if (min_num_fillin < 5) {
                            isIndependent = true;
                            break;
                        }
                    }
                }
                isIndependent = true;
            }
            if (min_num_fillin < 5) {
                break;
            }
        }

        if (isIndependent == false) {
            rank_Af = cnt;
            printf("rank_Af = %d\n", rank_Af);
            return;
        }

        setPivot(pivotRow, pivotCol, pivotVal, cnt);
        isPivoted[pivotCol] = true;

    }

#if LU_TRACE
    rMessage("");
#endif

    rank_Af = size_Kf;
    printf("rank_Af = %d, full-rank\n", rank_Af);
    return;
}


/* ============================================================
   getP
   ============================================================ */
int* LUFactor::getP()
{
    int* ret;
    NewArray(ret, int, rank_Af);
    for (int i = 0; i < rank_Af; i++) {
        ret[i] = P[i];
    }

    return ret;
}


/* ============================================================
   getQ
   ============================================================ */
int* LUFactor::getQ()
{
    int* ret;
    NewArray(ret, int, rank_Af);
    for (int i = 0; i < rank_Af; i++) {
        ret[i] = Q[i];
    }

    return ret;
}


/* ============================================================
   setRho
   ============================================================ */
void LUFactor::setRho(double new_rho)
{
    if( new_rho < 0 || new_rho > 1 ){
        rMessage( "LUInfo:setRho(): rho must be 0 <= rho <= 1" );
        return;
    }
    rho = new_rho;
    return;
}


/* ============================================================
   Destructor
   ============================================================ */
LUFactor::~LUFactor()
{
    if (pivot_Af != NULL) {
        delete pivot_Af;
        pivot_Af = NULL;
    }

    DeleteArray(P);
    DeleteArray(Q);
    for (int i = 0; i < size_m; i++) {
        DeleteArray(pos_Annz[i]);
    }
    DeleteArray(pos_Annz);
    DeleteArray(num_Annz);

    if (LInvP != NULL) {
        delete LInvP;
        LInvP = NULL;
    }

    DeleteArray(colPosition);
    return;
}
