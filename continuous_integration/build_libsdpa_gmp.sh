#!/bin/bash

# if [[ "$RUNNER_OS" == "RedHat" ]]; then
#     yum update -y
#     yum install -y epel-release
#     yum install -y glibc-static
# fi

if [[ "$RUNNER_OS" == "Windows" ]]; then
    cd /d
else
    cd $GITHUB_WORKSPACE
fi

# At every release of sdpa-python, we should attempt to link against latest
# GMP package (even if there is no new `sdpa-multiprecision` release)
# Please check the bottom of https://gmplib.org for latest GMP available.
# or http://ftp.de.debian.org/debian/pool/main/g/gmp
# NOTE: Debian does not have `docs` so archive will be slightly smaller
# we would however need to create a placeholder Makefile.in for it
curl -L -O http://ftp.de.debian.org/debian/pool/main/g/gmp/gmp_6.3.0+dfsg.orig.tar.xz
tar xf gmp_6.3.0+dfsg.orig.tar.xz
cd gmp-6.3.0+dfsg
mkdir doc
echo "all:" > doc/Makefile.in
echo "check:" >> doc/Makefile.in
echo "install:" >> doc/Makefile.in

if [[ "$RUNNER_OS" == "Windows" ]]; then
    # on Windows, longer paths cause problems during build
    # we therefore build GMP in /d
    ./configure --prefix=/d/gmp-6.3.0+dfsg --enable-cxx
    make
    make check
    make install
else
    if [[ "$RUNNER_OS" == "macOS" ]] && [[ "$RUNNER_ARCH" == "arm64" ]]; then
        ./configure --prefix=$GITHUB_WORKSPACE/gmp-6.3.0+dfsg --enable-cxx CFLAGS="-arch arm64" CXXFLAGS="-arch arm64" LDFLAGS="-arch arm64" --host=arm64-apple-darwin # arm64-apple-darwin23.2.0
    else # for Linux and macOS on x86_64
        ./configure --prefix=$GITHUB_WORKSPACE/gmp-6.3.0+dfsg --enable-cxx
    fi
    make
    make check
    make install
fi

ls -al lib # shows compiled libgmp
cd $GITHUB_WORKSPACE

if [[ "$RUNNER_OS" == "Windows" ]]; then
    cd sdpa-multiprecision
    ./configure --with-system-spooles --with-spooles-includedir=/d/msys64/mingw64/include/spooles CFLAGS="-DNDEBUG" CXXFLAGS="-DNDEBUG"
elif [[ "$RUNNER_OS" == "macOS" ]]; then
    cd sdpa-multiprecision
    if [[ "$RUNNER_ARCH" == "arm64" ]]; then
        sed -i.bak "s/-funroll-all-loops/-arch arm64 -funroll-all-loops/g" spooles/patches/patch-Make.inc
        ./configure --with-gmp-includedir=$GITHUB_WORKSPACE/gmp-6.3.0+dfsg/include --with-gmp-libdir=$GITHUB_WORKSPACE/gmp-6.3.0+dfsg/lib CFLAGS="-DNDEBUG -arch arm64" CXXFLAGS="-DNDEBUG -arch arm64" --host=arm64-apple-darwin
    else
        ./configure --with-gmp-includedir=$GITHUB_WORKSPACE/gmp-6.3.0+dfsg/include --with-gmp-libdir=$GITHUB_WORKSPACE/gmp-6.3.0+dfsg/lib CFLAGS="-DNDEBUG" CXXFLAGS="-DNDEBUG"
    fi
else
    cd sdpa-multiprecision
    ./configure --with-gmp-includedir=$GITHUB_WORKSPACE/gmp-6.3.0+dfsg/include --with-gmp-libdir=$GITHUB_WORKSPACE/gmp-6.3.0+dfsg/lib CFLAGS="-DNDEBUG" CXXFLAGS="-DNDEBUG"
fi

make
./sdpa_gmp example1.dat example1.out

cd $GITHUB_WORKSPACE

# Modify setupcfg.py file
if [[ "$RUNNER_OS" == "Windows" ]]; then
    original_value="${GITHUB_WORKSPACE}"
    new_value="${original_value//\\//}" # replace \ with /
    sed -i.bak 's@SDPA_DIR =.*@SDPA_DIR="'$new_value'\/sdpa-multiprecision"@g' sdpa-python/setupcfg.py
    sed -i.bak "s/MINGW_LIBS =.*/MINGW_LIBS=os.path.join('D:\/','msys64','mingw64','lib')/g" sdpa-python/setupcfg.py
    sed -i.bak "s/SPOOLES_INCLUDE =.*/SPOOLES_INCLUDE=os.path.join('D:\/','msys64','mingw64','include','spooles')/g" sdpa-python/setupcfg.py
    sed -i.bak "s/SPOOLES_DIR =.*/SPOOLES_DIR=os.path.join('D:\/','msys64','mingw64','lib')/g" sdpa-python/setupcfg.py
    sed -i.bak "s/GMP_DIR =.*/GMP_DIR=os.path.join('D:\/','gmp-6.3.0+dfsg')/g" sdpa-python/setupcfg.py
    sed -i.bak "s/USEGMP =.*/USEGMP = True/g" sdpa-python/setupcfg.py
    echo "[build]" > sdpa-python/setup.cfg
    echo "compiler=mingw32" >> sdpa-python/setup.cfg
elif [[ "$RUNNER_OS" == "macOS" ]]; then
    sed -i.bak 's@SDPA_DIR =.*@SDPA_DIR="'"$GITHUB_WORKSPACE"'/sdpa-multiprecision"@g' sdpa-python/setupcfg.py
    sed -i.bak 's@GMP_DIR =.*@GMP_DIR="'"$GITHUB_WORKSPACE"'/gmp-6.3.0+dfsg"@g' sdpa-python/setupcfg.py
    sed -i.bak "s/USEGMP =.*/USEGMP = True/g" sdpa-python/setupcfg.py
else
    sed -i.bak 's@SDPA_DIR =.*@SDPA_DIR="'"$GITHUB_WORKSPACE"'/sdpa-multiprecision"@g' sdpa-python/setupcfg.py
    sed -i.bak 's@GMP_DIR =.*@GMP_DIR="'"$GITHUB_WORKSPACE"'/gmp-6.3.0+dfsg"@g' sdpa-python/setupcfg.py
    sed -i.bak "s/USEGMP =.*/USEGMP = True/g" sdpa-python/setupcfg.py
fi
