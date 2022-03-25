/* ********************************************************************************
   fvelim_LUFactor.h
   A module of fvelim
   C interface to call free variable elimination

   December 2010, Kenta KATO
   ******************************************************************************** */

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
