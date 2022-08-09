/*
C Extension code to call SPOOLES
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

September 2010: Originally written by Kenta Kato
*/

extern "C" {
// Spooles include
#include <InpMtx.h>
#include <IVL.h>
#include <Graph.h>
#include <ETree.h>
#include <misc.h>
#include <IV.h>
#include <SymbFac.h>
};
#include <iostream>
#include <cstdio>
#include "spcoloextmodule.h"
#include "spcolo_util.h"

#define MX_DEBUG 0
#define SPARSE_DEBUG 1

#if SPARSE_DEBUG
#define AGGREGATE_THRESHOLD 2.0
#define EXTEND_THRESHOLD 2.0
#else
#define AGGREGATE_THRESHOLD 0.25
#define EXTEND_THRESHOLD 0.4
#endif



int _countNonZero(int m, IVL* symbfacIVL);

/* ============================================================
   sclo_ordering_mmd
   Get the permutation ordering by minimum degree ordering algorithm
   ============================================================ */
int* spcolo_ordering_mmd(int* rowind, int* colptr, int size_n)
{
    /* Initialize Spooles objects */
    IVL* adjIVL = IVL_new();
    Graph* graph = Graph_new();

    /* --------------------------------------------------
       Construct adjacency list of SPOOLES
       -------------------------------------------------- */
    IVL_init1(adjIVL, IVL_CHUNKED, size_n);

    int size_vec;
    int* vec_ind;
    NewArray(vec_ind, int, size_n);
    /*
    int* vec_ind = new int[size_n];
    */
    for (int col = 0; col < size_n; col++) {
        size_vec = 0;
        for (int row = colptr[col]; row < colptr[col + 1]; row++) {
            vec_ind[size_vec] = rowind[row];
            size_vec++;
        }
        IVL_setList(adjIVL, col, size_vec, vec_ind);
    }

    /* --------------------------------------------------
       Construct graph of SPOOLES
       -------------------------------------------------- */
    Graph_init2(graph, 0, size_n, 0, IVL_tsize(adjIVL), size_n, IVL_tsize(adjIVL),
                adjIVL, NULL, NULL);

    DeleteArray(vec_ind);

    if (IVL_tsize(adjIVL) > AGGREGATE_THRESHOLD * size_n * size_n) {
        /* Not need to decompose matrix */
        Graph_free(graph);
        return NULL;
    }

    /* --------------------------------------------------
       Spooles MMD
       -------------------------------------------------- */
    ETree* etree;
    int seed = 0;
    int msglvl = 0;
    FILE* fp = NULL;
    IV* newToOldIV_MMD;
    IVL* symbfacIVL_MMD;
    etree = orderViaMMD(graph, seed, msglvl, fp);
    newToOldIV_MMD = ETree_newToOldVtxPerm(etree);
    symbfacIVL_MMD = SymbFac_initFromGraph(etree, graph);

    int nonzeros = _countNonZero(size_n, symbfacIVL_MMD) * 2 - size_n;

    ETree_free(etree);
    Graph_free(graph);

    if (nonzeros > EXTEND_THRESHOLD * size_n * size_n) {
        /* Not need to decompose matrix */
        return NULL;
    }

    /* --------------------------------------------------
       Get ordering list and copy to PyObject
       -------------------------------------------------- */
    return IV_entries(newToOldIV_MMD);
}


/* ============================================================
   Count number of nonzero of IVL
   ============================================================ */
int _countNonZero(int m, IVL* symbfacIVL)
{
    int nonzeros = 0;
    bool bnode[m];

    // count nonzero element
    for (int i = 0; i < m; i++) {
        bnode[i] = false;
    }

    int nClique = IVL_nlist(symbfacIVL);
    int psize;
    int* pivec;
    for (int l = nClique - 1; l >= 0; l--) {
        IVL_listAndSize(symbfacIVL, l, &psize, &pivec);
        for (int i = 0; i < psize; i++) {
            int ii = pivec[i];
            if (bnode[ii] == false) {
                nonzeros += psize - i;
                bnode[ii] = true;
            }
        }
    }

    return nonzeros;
}
