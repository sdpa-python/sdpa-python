import os
import platform

"""
If you intend to use `sdpa-multiprecision` backend, set this to True
"""
USEGMP = False

"""
>> REQUIRED ALWAYS

Specify the package name of the backend you intend to use.
The regular backend can be downloaded from https://sdpa.sourceforge.net
and currently the latest is `sdpa-7.3.16`

The latest version of the multiprecision backend can be cloned from
https://github.com/sdpa-python/sdpa-multiprecision
"""
sdpa_package_name = 'sdpa-7.3.16'
# OR sdpa_package_name = 'sdpa-multiprecision'

"""
>> REQUIRED ALWAYS

Provide the location SDPA libraries and headers

(Windows template)
SDPA_DIR = os.path.join("C:\\", "path", "to", sdpa_package_name)

(Linux, macOS template)
SDPA_DIR = f'/path/to/{sdpa_package_name}'
"""

SDPA_DIR = ""


"""
(AUTOSET) MUMPS is only used by the regular backend. SDPA buildsystem
is configured to download and build MUMPS as part of the SDPA package. 
"""

MUMPS_DIR =  os.path.join(SDPA_DIR, "mumps", "build")

"""
>> REQUIRED ALWAYS

The `spcolo` module of this Python frontend requires SPOOLES.
If using the multiprecision backend, the backend also uses SPOOLES.
The multiprecision backend contains SPOOLES. If using `sdpa-multiprecision`
on Linux or macOS, leave this unchanged. However, if you are using
the regular backend, you can download SPOOLES from
http://www.netlib.org/linalg/spooles/spooles.2.2.html

(Windows template)
SPOOLES_DIR = os.path.join("C:\\", "path", "to", "spooles")

(Linux, macOS template)
SPOOLES_DIR = '/path/to/spooles/'

On Windows, if you use SPOOLES provided by MSYS2, you can set 
SPOOLES_INCLUDE = os.path.join("C:\\", "msys64", "mingw64", "include", "spooles")
SPOOLES_DIR need not be changed as `libspooles.a` is in `MINGW_LIBS` path and
is accessible to setup.py
"""

SPOOLES_INCLUDE = os.path.join(SDPA_DIR, "spooles", "build")
SPOOLES_DIR =  os.path.join(SDPA_DIR, "spooles", "build")

"""
>> REQUIRED ONLY: if building on Windows

This probably does not require change if you installed MSYS2 in
the default location suggested by the installer.
"""
MINGW_LIBS =  os.path.join("C:\\", "msys64", "mingw64", "lib")

"""
>> REQUIRED ONLY: for the regular backend

On Linux or Windows , we can use either Reference BLAS
BLAS_LAPACK_LIBS = ['blas', 'lapack']

Or we can use OpenBLAS (recommended)
BLAS_LAPACK_LIBS = ['openblas']

If using OpenBLAS on Windows, you might need to link against `libgomp.a`
BLAS_LAPACK_LIBS = ['openblas', 'gomp']

On macOS, we use the BLAS provided by Accelerate
BLAS_LAPACK_LIBS = []
"""

BLAS_LAPACK_LIBS = []

"""
>> REQUIRED ONLY: for regular backend on macOS

`gfortran` builds MUMPS (used by regular backend) and the location of `libgfortran.a`
is in a known path on both Linux and MinGW (Windows).

On macOS, you need to specify the location of `libgfortran.a`

You can find the library using `find / -name 'libgfortran.a'`

If you install using HomeBrew, currently it looks like
GFORTRAN_LIBS = '/usr/local/Cellar/gcc/12.2.0/lib/gcc/current'
"""
GFORTRAN_LIBS = '/usr/local/Cellar/gcc/12.2.0/lib/gcc/current'

"""
>> REQUIRED ONLY: for multiprecision backend

Provide the location the GNU multiprecision library

(Windows template)
GMP_DIR = os.path.join("C:\\", "path", "to", f"gmp-{version}")

(Linux, macOS template)
GMP_DIR = f'/path/to/gmp-{version}'

Current latest version is 6.2.1
"""

version = '6.2.1'
GMP_DIR = ""
