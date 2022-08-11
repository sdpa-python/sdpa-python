import os
import platform

# SPOOLES can be downloaded from http://www.netlib.org/linalg/spooles/spooles.2.2.html

version = '7.3.16'

# This will need to be set if building on Linux or macOS
SDPA_MAKEINC  = f'/path/to/sdpa-{version}/etc/make.inc'
SPOOLES_INCLUDE = '/path/to/spooles/' # Linux, macOS template

# This will need to be set ONLY if building on macOS
GFORTRAN_LIBS ='/usr/local/gfortran/lib'

# These will need to be set ONLY if building on Windows
SDPA_DIR      =  os.path.join("C:\\", "path", "to", f"sdpa-{version}")
# MINGW_LIBS    =  os.path.join("C:\\", "msys64", "mingw64", "lib")
# SPOOLES_INCLUDE =  os.path.join("C:\\", "path", "to", "spooles") # Windows template

# (Optional) paths for Lapack and Blas
LAPACK_DIR = '/optional/path/to/lapack/' # library may already be installed system wide
LAPACK_NAME = 'lapack'
BLAS_DIR = '/optional/path/to/blas' # library may already be installed system wide
BLAS_NAME = 'blas'

# No changes needed below

SPOOLES_DIR = SPOOLES_INCLUDE
SPOOLES_NAME = 'spooles'

if platform.system()=='Windows':
  MUMPS_DIR   =  os.path.join(SDPA_DIR, "mumps", "build")
