#!/bin/bash

if [[ "$RUNNER_OS" == "Linux" ]]; then
    # when this script runs directly inside a Ubuntu runner
    sudo apt-get update -y
    sudo apt-get install -y libopenblas-dev
elif [[ "$RUNNER_OS" == "manylinux" ]]; then
    # when this script runs inside a manylinux container
    yum update -y
    yum install -y epel-release
    yum install -y blas-devel lapack-devel
    yum install -y openblas-devel
fi

# Download and patch libsdpa
curl -L -O https://downloads.sourceforge.net/project/sdpa/sdpa/sdpa_7.3.17.tar.gz
tar -zxf sdpa_7.3.17.tar.gz
cd sdpa-7.3.17
# At every release of sdpa-python, we should attempt to link against latest
# MUMPS package (even if there is no new `sdpa` release at sdpa.sourceforge.net)
# Please check the bottom of this webpage for latest MUMPS available:
# http://ftp.de.debian.org/debian/pool/main/m/mumps/
sed -i.bak 's/MUMPS_VER =.*/MUMPS_VER = 5.6.2/' mumps/Makefile
# `wget` may not be available but `curl` almost always is
sed -i.bak 's/wget/curl -L -O/' mumps/Makefile

if [[ "$RUNNER_OS" == "macOS" ]]; then
    ./configure CFLAGS="-DNDEBUG" CXXFLAGS="-DNDEBUG" FCFLAGS="-DNDEBUG"
else
    ./configure --with-blas="-lopenblas" --with-lapack="-lopenblas" CFLAGS="-DNDEBUG" CXXFLAGS="-DNDEBUG" FCFLAGS="-DNDEBUG"
fi

# Build SDPA and MUMPS and `cd` back into main directory
make
cd $GITHUB_WORKSPACE

# Build SPOOLES or provide location to system spooles
if [[ "$RUNNER_OS" != "Windows" ]]; then
    # Download and build SPOOLES
    curl -L -O http://ftp.de.debian.org/debian/pool/main/s/spooles/spooles_2.2.orig.tar.gz
    mkdir spooles
    tar -zxf spooles_2.2.orig.tar.gz -C spooles
    cd spooles
    sed -i.bak 's/CC =.*//' Make.inc
    if [[ "$RUNNER_OS" == "Linux" ]] || [[ "$RUNNER_OS" == "manylinux" ]]; then
        sed -i.bak 's/^  CFLAGS =.*/& -fPIC/' Make.inc
        # this patch is critical only on Ubuntu
        curl -O https://raw.githubusercontent.com/sdpa-python/sdpa-multiprecision/main/spooles/patches/patch-timings.h
        patch -p0 < patch-timings.h
    fi
    make lib
    mv spooles.a libspooles.a
    cd $GITHUB_WORKSPACE
fi

# Modify setupcfg.py file
if [[ "$RUNNER_OS" == "Windows" ]]; then
    original_value="${GITHUB_WORKSPACE}"
    new_value="${original_value//\\//}" # replace \ with /
    sed -i.bak 's@SDPA_DIR =.*@SDPA_DIR="'$new_value'\/sdpa-7.3.17"@g' sdpa-python/setupcfg.py
    sed -i.bak "s/MINGW_LIBS =.*/MINGW_LIBS=os.path.join('D:\/','msys64','mingw64','lib')/g" sdpa-python/setupcfg.py
    sed -i.bak "s/SPOOLES_INCLUDE =.*/SPOOLES_INCLUDE=os.path.join('D:\/','msys64','mingw64','include','spooles')/g" sdpa-python/setupcfg.py
    sed -i.bak "s/SPOOLES_DIR =.*/SPOOLES_DIR=os.path.join('D:\/','msys64','mingw64','lib')/g" sdpa-python/setupcfg.py
    sed -i.bak "s/BLAS_LAPACK_LIBS =.*/BLAS_LAPACK_LIBS = ['openblas', 'gomp']/g" sdpa-python/setupcfg.py
    echo "[build]" > sdpa-python/setup.cfg
    echo "compiler=mingw32" >> sdpa-python/setup.cfg
elif [[ "$RUNNER_OS" == "macOS" ]]; then
    sed -i.bak 's@SDPA_DIR =.*@SDPA_DIR="'"$GITHUB_WORKSPACE"'/sdpa-7.3.17"@g' sdpa-python/setupcfg.py
    sed -i.bak 's@SPOOLES_DIR =.*@SPOOLES_DIR="'"$GITHUB_WORKSPACE"'/spooles"@g' sdpa-python/setupcfg.py
    sed -i.bak 's@SPOOLES_INCLUDE =.*@SPOOLES_INCLUDE="'"$GITHUB_WORKSPACE"'/spooles"@g' sdpa-python/setupcfg.py
    # check if /usr/local/opt/gcc/lib/gcc/13 is more generic
    sed -i.bak "s/GFORTRAN_LIBS =.*/GFORTRAN_LIBS='\/usr\/local\/Cellar\/gcc\/13.2.0\/lib\/gcc\/current'/g" sdpa-python/setupcfg.py
else
    sed -i.bak 's@SDPA_DIR =.*@SDPA_DIR="'"$GITHUB_WORKSPACE"'/sdpa-7.3.17"@g' sdpa-python/setupcfg.py
    sed -i.bak 's@SPOOLES_DIR =.*@SPOOLES_DIR="'"$GITHUB_WORKSPACE"'/spooles"@g' sdpa-python/setupcfg.py
    sed -i.bak 's@SPOOLES_INCLUDE =.*@SPOOLES_INCLUDE="'"$GITHUB_WORKSPACE"'/spooles"@g' sdpa-python/setupcfg.py
    sed -i.bak "s/BLAS_LAPACK_LIBS =.*/BLAS_LAPACK_LIBS = ['openblas']/g" sdpa-python/setupcfg.py

    # CentOS 7 aarch does not have libquadmath and does not need it
    # libquadmath is an indirect dependency, due to libgfortran
    # for x86_64, repairwheel will automatically static link quadmath
    sed -i.bak "s/, 'quadmath'//" sdpa-python/setup.py
fi
