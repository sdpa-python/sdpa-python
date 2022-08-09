#!/usr/bin/env python
"""
Class definition of MatData (used as interface between SciPy and SDPA C API)
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
"""

from scipy import sparse, array

class MatData(object):
    """Description of a sparse matrix

    Attributes:
      size_row: Size of row
      size_col: Size of column
      values: List of Values
      rowind: List of row indices
      colptr: List of start index of each column
    """
    def __init__(self, mat=None,
                 values=None, rowind=None, colptr=None, size=None):
        """Constructor for sparse matrix

        Args:
          mat: input sparse matrix
          values: List of Values
          rowind: List of row indices
          colptr: List of start index of each column
          size: Tuple (size_row, size_col)
        """
        if sparse.issparse(mat):
            mat2 = mat.tocsc() if not sparse.isspmatrix_csc(mat) else mat
            self.size_row, self.size_col = mat2.shape
            self.values = list(mat2.data)
            self.rowind = list(mat2.indices)
            self.colptr = list(mat2.indptr)
        elif (isinstance(values, list) and isinstance(rowind, list) and
              isinstance(colptr, list) and isinstance(size, tuple)):
            if (len(values) != len(rowind) or len(colptr) != size[1] + 1):
                raise ValueError('size of val, row, col list must be same')
            self.size_row, self.size_col = size
            self.values, self.rowind, self.colptr = values, rowind, colptr
        else:
            raise TypeError('Input arg must be sparse matrix or lists.')

    def __len__(self):
        """Return num of nonzero elements in matrix"""
        return len(self.values)

    def size(self):
        """Return the size of matrix"""
        return (self.size_row, self.size_col)

    def tomatrix(self):
        """Get sparse matrix from self

        Returns:
          Sparse Matrix in lil_matrix format
        """
        return sparse.csc_matrix((array(self.values),
                                  array(self.rowind),
                                  array(self.colptr)),
                                 shape=(self.size_row, self.size_col))
