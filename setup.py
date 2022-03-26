from setuptools import setup, Extension

from setupcfg import *

pathdict = dict((s[0].strip(), s[1].strip()) for s in
            [line.split('=') for line in open(SDPA_MAKEINC).readlines()])

SDPA_DIR = pathdict['SDPA_DIR']
SDPA_LIB = SDPA_DIR # + '/lib'
SDPA_INCLUDE = SDPA_DIR + '/include'

MUMPS_DIR = pathdict['MUMPS_DIR']
MUMPS_LIB = MUMPS_DIR + '/lib'
MUMPS_LIBSEQ = MUMPS_DIR + '/libseq'
MUMPS_INCLUDE = MUMPS_DIR + '/include'

if not LAPACK_NAME:
    LAPACK_NAME = 'lapack'
if not BLAS_NAME:
    BLAS_NAME = 'blas'
if not SPOOLES_INCLUDE:
    SPOOLES_INCLUDE = '/usr/include/spooles/'
if not SPOOLES_NAME:
    SPOOLES_NAME = 'spooles'


ext_sdpacall = Extension(
    'sdpap.sdpacall.sdpa',
    [
        'sdpap/sdpacall/cmodule/sdpamodule.cpp'
    ],
    include_dirs=[SDPA_DIR, MUMPS_INCLUDE],
    library_dirs=[SDPA_LIB, MUMPS_LIB, MUMPS_LIBSEQ], #LAPACK_DIR, BLAS_DIR],
    libraries=['sdpa', 'dmumps', 'mumps_common', 'pord', 'mpiseq', LAPACK_NAME, BLAS_NAME]
)

ext_fvelim = Extension(
    'sdpap.fvelim.fvelimext',
    [
        'sdpap/fvelim/cmodule/fvelimextmodule.cpp',
        'sdpap/fvelim/cmodule/fvelim_LUFactor.cpp',
        'sdpap/fvelim/cmodule/fvelim_SparseMatrix.cpp'
    ],
    include_dirs=['sdpap/fvelim/cmodule'],
)

ext_spcolo = Extension(
    'sdpap.spcolo.spcoloext',
    [
        'sdpap/spcolo/cmodule/spcoloextmodule.cpp',
        'sdpap/spcolo/cmodule/spcolo_cholesky.cpp',
        'sdpap/spcolo/cmodule/spcolo_SparseMatrix.cpp',
        'sdpap/spcolo/cmodule/spcolo_ordering.cpp'
    ],
    include_dirs=['sdpap/spcolo/cmodule', SPOOLES_INCLUDE],
    library_dirs=[],#[SPOOLES_DIR],
    libraries=[SPOOLES_NAME],
)

setup(
    name='sdpa-python',
    version='0.0.1',
    description='Python 3 Wrapper for SDPA (SemiDefinite Programming Algorithm)',
    long_description=open('README.md').read(),
    author='Usama Muneeb',
    packages=[
        'sdpap',
        'sdpap.sdpacall',
        'sdpap.fvelim',
        'sdpap.spcolo'
    ],
    ext_modules=[
        ext_sdpacall,
        ext_fvelim,
        ext_spcolo
    ],
    python_requires=">=3.6"
)
