#!/usr/bin/env python
"""
Class definition of `CliqueSet`
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
"""

__all__ = ['CliqueSet']

class CliqueSet(object):
    """Maximal clique set of Aggregate sparsity pattern

    Attributes:
      cliques: A list of clique. Each clique is described as list of integer.
      map_convIndex: A dictionary of mapping list from index of ASP (i,j)
        to converted index.
        i.e. map_convIndex[(i,j)] = converted_index
      map_origIndex: A dictionary of mapping list from converted index
        to original index.
        i.e. map_index[converted_index] = original_index
    """
    def __init__(self, size_n):
        """Constructor"""
        self.cliques = []
        self.map_convIndex = {}
        self.map_origIndex = {}
        self.num_element = 0
        self.size_n = size_n

    def __len__(self):
        """Return number of cliques"""
        return len(self.cliques)

    def append_clique(self, clq):
        """Append one clique to cliques and make index mapping

        Args:
          clq: A list of integer.
        """
        self.cliques.append(clq)
        for i in clq:
            for j in clq:
                if i >= j:
                    origIndex1 = i * self.size_n + j
                    origIndex2 = j * self.size_n + i
                    if origIndex1 not in self.map_convIndex:
                        self.map_origIndex[self.num_element] = \
                                                        (origIndex1, origIndex2)
                        self.map_convIndex[origIndex1] = self.num_element
                        self.map_convIndex[origIndex2] = self.num_element
                        self.num_element += 1
        return
