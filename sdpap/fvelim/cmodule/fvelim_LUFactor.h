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

#ifndef __FVELIM_LUFACTOR_H__
#define __FVELIM_LUFACTOR_H__

#include "fvelim_SparseMatrix.h"
#include "fvelim_util.h"

using namespace std;

/* ============================================================
   LUFactor
   LU decomposed matrix information
   ============================================================ */
class LUFactor
{
 private:

    /* Pivotted free variable matrix */
    SparseMatrix* pivot_Af;

    /* Size of Af's col */
    const int size_Kf;

    /* Size of (Al, Aq, As)'s col */
    const int size_Alqs;

    /* Size of Af's row or matrix T */
    const int size_m;

    /* Matrix rank of Af */
    int rank_Af;

    /* Permutation from left */
    int* P;

    /* Permutation from right */
    int* Q;

    /* Position of nonzero in pivotted (Al, Aq, As) */
    int** pos_Annz;

    /* Number of nonzero in each pivotted row of (Al, Aq, As) */
    int* num_Annz;

    /* Matrix L^(-1) * P */
    SparseMatrix* LInvP;

    /* Parameter for selecting Marcowitz candidate */
    double rho;

    /* Zero point for numeric error */
    double zeroValue;

    /*
     * To make P, Q
     */

    /* Col index position after pivot */
    int* colPosition;

    /* --------------------------------------------------
       Update pivot_Af or LInvP
       -------------------------------------------------- */
    void _updateMatrix(SparseMatrix* target, SparseMatrix* T,
                       int pivotRow, int pivotCol, int index);

 public:
    /* --------------------------------------------------
       Constructor
       -------------------------------------------------- */
    LUFactor(SparseMatrix* A_f, SparseMatrix* A_lqs, double rho, double zeroPoint);

    /* --------------------------------------------------
       setPivot
       -------------------------------------------------- */
    void setPivot(int pivotRow, int pivotCol, double pivotVal, int index);

    /* --------------------------------------------------
       getU
       -------------------------------------------------- */
    SparseMatrix* getU();

    /* --------------------------------------------------
       getV
       not full rank case
       -------------------------------------------------- */
    SparseMatrix* getV();

    /* --------------------------------------------------
       Get L^(-1) * P
       -------------------------------------------------- */
    SparseMatrix* getLInvP();

    /* --------------------------------------------------
       decompose
       LU Decomposition
       -------------------------------------------------- */
    void decompose();

    /* --------------------------------------------------
       getP
       Get permutation from left
       -------------------------------------------------- */
    int* getP();

    /* --------------------------------------------------
       getQ
       Get permutation from right
       -------------------------------------------------- */
    int* getQ();

    /* --------------------------------------------------
       setRho
       Set parameter rho
       -------------------------------------------------- */
    void setRho(double new_rho);

    /* --------------------------------------------------
       getRank
       Get matrix rank Af
       -------------------------------------------------- */
    int getRank(){ return rank_Af; }

    /* --------------------------------------------------
       getRho
       Get parameter rho
       -------------------------------------------------- */
    double getRho(){ return rho; }

    /* --------------------------------------------------
       Destructor
       -------------------------------------------------- */
    ~LUFactor();
};

#endif
