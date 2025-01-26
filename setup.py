from setuptools import setup, Extension
import platform

from setupcfg import *

def pjoin(parent, children):
    return list(map(lambda x : os.path.join(parent, x), children))

if USEGMP:
    include_dirs = pjoin(SDPA_DIR, ['.', 'mpack']) + [SPOOLES_INCLUDE] + pjoin(GMP_DIR, ['include'])
    libraries = ['sdpa_gmp', 'spooles']
    library_dirs = [SDPA_DIR, SPOOLES_DIR]
    if False: # platform.system()=='Darwin':
        extra_objects = pjoin(GMP_DIR, ['lib/libgmpxx.a'])
    else:
        libraries += ['gmp', 'gmpxx']
        library_dirs += pjoin(GMP_DIR, ['lib'])
        extra_objects = []
else:
    libraries = ['sdpa', 'dmumps', 'mumps_common', 'pord', 'mpiseq'] + BLAS_LAPACK_LIBS
    library_dirs = [SDPA_DIR] + pjoin(MUMPS_DIR, ['lib', 'libseq'])
    include_dirs = [SDPA_DIR] + pjoin(MUMPS_DIR, ['include'])
    if platform.system()=='Darwin':
        extra_objects = pjoin(GFORTRAN_LIBS, ['libgfortran.a', 'libquadmath.a'])
    else:
        libraries += ['gfortran', 'quadmath']
        extra_objects = []

if platform.system()=='Windows':
    import distutils.cygwinccompiler
    distutils.cygwinccompiler.get_msvcr = lambda: []
    library_dirs += [MINGW_LIBS]

ext_sdpacall = Extension(
    'sdpap.sdpacall.sdpa',
    [
        'sdpap/sdpacall/cmodule/sdpamodule.cpp'
    ],
    include_dirs=include_dirs,
    library_dirs=library_dirs,
    libraries=libraries,
    extra_objects=extra_objects,
    extra_compile_args=['-DUSEGMP=1'] if USEGMP else ['-DUSEGMP=0'],
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
    library_dirs=[SPOOLES_DIR],
    libraries=['spooles'],
    extra_link_args=['-static'] if platform.system()=='Windows' else []
)

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='sdpa-multiprecision' if USEGMP else 'sdpa-python',
    version='0.2.2',
    maintainer='Usama Muneeb',
    maintainer_email='umunee2@uic.edu',
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
    python_requires=">=3.7",

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
