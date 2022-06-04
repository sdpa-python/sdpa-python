import os

# This will need to be set if building on Linux or macOS
SDPA_MAKEINC = '/path/to/sdpa/etc/make.inc'

# This will need to be set if building on Windows
SDPA_DIR    =  os.path.join("C:\\", "path", "to", "sdpa-7.3.8")


MUMPS_DIR   =  os.path.join(SDPA_DIR, "mumps", "build")


# (Optional) paths for Lapack and Blas
LAPACK_DIR = '/optional/path/to/lapack/' # library may already be installed system wide
LAPACK_NAME = 'lapack'
BLAS_DIR = '/optional/path/to/blas' # library may already be installed system wide
BLAS_NAME = 'blas'

# Please download SPOOLES from http://www.netlib.org/linalg/spooles/spooles.2.2.html
# While library may be installed system wide, headers must be provided

SPOOLES_INCLUDE = '/path/to/spooles/' # Linux, macOS template
# SPOOLES_INCLUDE =  os.path.join("C:\\", "path", "to", "spooles") # Windows template

SPOOLES_DIR = SPOOLES_INCLUDE
SPOOLES_NAME = 'spooles'
