from setuptools import setup, Extension
import platform

from setupcfg import *

if platform.system()=='Windows':
    import distutils.cygwinccompiler
    distutils.cygwinccompiler.get_msvcr = lambda: []
else:
    pathdict = dict((s[0].strip(), s[1].strip()) for s in
                [line.split('=') for line in open(SDPA_MAKEINC).readlines()])

    SDPA_DIR = pathdict['SDPA_DIR']
    MUMPS_DIR = pathdict['MUMPS_DIR']


SDPA_LIB = SDPA_DIR # + '/lib'
SDPA_INCLUDE = SDPA_DIR + '/include'

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

EXTRA_LIBS = []
if platform.system()=='Windows':
    EXTRA_LIBS = [MINGW_LIBS]
# elif platform.system()=='Darwin':
#     EXTRA_LIBS = [GFORTRAN_LIBS]

STATIC_LIBS = ['gfortran', 'quadmath']
if platform.system()=='Darwin':
    STATIC_LIBS = []
    GFORTRAN_LIBS_STATIC = list(map(lambda x : os.path.join(GFORTRAN_LIBS, x), ['libgfortran.a', 'libquadmath.a']))

ext_sdpacall = Extension(
    'sdpap.sdpacall.sdpa',
    [
        'sdpap/sdpacall/cmodule/sdpamodule.cpp'
    ],
    include_dirs=[SDPA_DIR, MUMPS_INCLUDE],
    library_dirs=[SDPA_LIB, MUMPS_LIB, MUMPS_LIBSEQ] + EXTRA_LIBS,
    libraries=['sdpa', 'dmumps', 'mumps_common', 'pord', 'mpiseq', LAPACK_NAME, BLAS_NAME] + STATIC_LIBS,
    extra_objects=GFORTRAN_LIBS_STATIC if platform.system()=='Darwin' else [],
    extra_link_args=['-static'] if platform.system()=='Windows' else []
)

ext_fvelim = Extension(
    'sdpap.fvelim.fvelimext',
    [
        'sdpap/fvelim/cmodule/fvelimextmodule.cpp',
        'sdpap/fvelim/cmodule/fvelim_LUFactor.cpp',
        'sdpap/fvelim/cmodule/fvelim_SparseMatrix.cpp'
    ],
    include_dirs=['sdpap/fvelim/cmodule'],
    extra_link_args=['-static'] if platform.system()=='Windows' else []
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
    library_dirs=[SPOOLES_DIR] if platform.system()=='Darwin' else [],
    libraries=[SPOOLES_NAME],
    extra_link_args=['-static'] if platform.system()=='Windows' else []
)

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='sdpa-python',
    version='0.1',
    maintainer='Usama Muneeb',
    url='https://sdpa-python.github.io',
    description='SDPA (SemiDefinite Programming Algorithm) for Python',
    long_description=long_description,
    long_description_content_type="text/markdown",

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
    python_requires=">=3.6",

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',

    ],
    license = 'GNU GPL version 2',
    install_requires=[
        "scipy >= 1.1.0"
    ],
)
