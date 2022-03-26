# setupcfg.py
# SDPAP configuration file
#
# Change the paths below depending on your settings.

# Path for make.inc of SDPA
SDPA_MAKEINC = '/path/to/sdpa/etc/make.inc'

# (Optional) paths for Lapack and Blas
LAPACK_DIR = '/optional/path/to/lapack/' # library may already be installed system wide
LAPACK_NAME = 'lapack'
BLAS_DIR = '/optional/path/to/blas' # library may already be installed system wide
BLAS_NAME = 'blas'

# Path for Spooles
SPOOLES_DIR = '/optional/path/to/spooles/' # library may already be installed system wide
SPOOLES_INCLUDE = '/path/to/spooles/' # headers must be provided
                                      # please download from http://www.netlib.org/linalg/spooles/spooles.2.2.html
SPOOLES_NAME = 'spooles'
