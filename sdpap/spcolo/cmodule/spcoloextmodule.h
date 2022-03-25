/* ************************************************************
   spcoloextmodule.h
   A module of spcolo

   October 2010, Kenta KATO
   ************************************************************ */

#ifndef __SPCOLOEXTMODULE_H__
#define __SPCOLOEXTMODULE_H__

#include "spcolo_SparseMatrix.h"

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define max(x, y) ((x) >= (y) ? (x) : (y))

using namespace std;

/* ============================================================
   Function declaration
   ============================================================ */

/* Get the permutation ordering by MMD algorithm */
int* spcolo_ordering_mmd(int* rowind, int* colptr, int size_n);

/* Cholesky factorization */
SparseMatrix* spcolo_cholesky(SparseMatrix* A);

#endif
