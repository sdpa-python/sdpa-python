/*
Basic required by SparseCoLO subprograms
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

October 2010: Originally written by Kenta Kato
*/

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
