"""
Copyright (C) 2023 Usama Muneeb

This file checks the consistency of the read/write methods that
import/export problem data from/to SDPA Sparse and CLP formatted files.

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
"""

import os
import sdpap

def check_clp_consistency(reference_filename, accuracy="%+8.16e"):
    reference_src = os.path.join("examples", reference_filename)
    output_dir = os.path.join(".", "tmp")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    output_src = os.path.join(output_dir, reference_filename)

    A, b, c, K, J = sdpap.readproblem(reference_src)
    sdpap.writeproblem(output_src, A, b, c, K, J, accuracy)
    A2, b2, c2, K2, J2 = sdpap.readproblem(output_src)

    diff_A = abs(A2-A).max()
    diff_B = abs(b2-b).max()
    diff_C = abs(c2-c).max()

    epsilon = 1e-30 # accounting for small numerical problems
    assert (diff_A < epsilon and diff_B < epsilon and diff_C < epsilon)


def check_sdpa_consistency(reference_filename, accuracy="%+8.16e"):
    reference_src = os.path.join("examples", reference_filename)
    output_dir = os.path.join(".", "tmp")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    output_src = os.path.join(output_dir, reference_filename)

    A, b, c, K, J = sdpap.fromsdpa(reference_src)
    sdpap.tosdpa(output_src, A, b, c, K, J, accuracy)
    A2, b2, c2, K2, J2 = sdpap.fromsdpa(output_src)

    diff_A = abs(A2-A).max()
    diff_B = abs(b2-b).max()
    diff_C = abs(c2-c).max()

    epsilon = 1e-30 # accounting for small numerical problems
    assert (diff_A < epsilon and diff_B < epsilon and diff_C < epsilon)


class TestReadWriteCLP():
    def test_clp_consistency_1(self):
        return check_clp_consistency('example1.clp')

    def test_clp_consistency_2(self):
        return check_clp_consistency('example2.clp', '%+1.6e')

    def test_clp_consistency_3(self):
        return check_clp_consistency('example3.clp')


class TestReadWriteSDPA():
    # this is provided with SDPA source code as well as user manuals
    def test_sdpa_consistency_1(self):
        return check_sdpa_consistency('sdpa_example.dat-s')

    # the following are from SDPLIB
    def test_sdpa_consistency_2(self):
        return check_sdpa_consistency('control2.dat-s')

    def test_sdpa_consistency_3(self):
        return check_sdpa_consistency('mcp500-4.dat-s')

    def test_sdpa_consistency_4(self):
        return check_sdpa_consistency('theta6.dat-s')

    # This contains some zero rows
    def test_sdpa_consistency_5(self):
        return check_sdpa_consistency('ss30.dat-s')

