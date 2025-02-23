#!/bin/bash

INSTALL_DIR=$(pwd)/external/install
mkdir -p $INSTALL_DIR

###############################################################################
#                                HEADER-ONLY
###############################################################################
#   eigen
#   spdlog

###############################################################################
#                                   SUNDIALS
###############################################################################
cd ./external/sundials
mkdir -p build
cd build

# BUILD_TESTING is only used when building examples is enabled
cmake_flags="-DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/sundials \
             -DBUILD_SHARED_LIBS=OFF \
             -DBUILD_TESTING=OFF \
             -DEXAMPLES_INSTALL=OFF \
             -DEXAMPLES_ENABLE_C=OFF \
             -DEXAMPLES_ENABLE_CXX=OFF"

if [[ "$(uname)" == "Linux" ]]; then
    cmake_flags="$cmake_flags -DCMAKE_C_FLAGS=-fPIC"
elif [[ "$(uname)" == "Darwin" ]]; then
    cmake_flags="$cmake_flags -DCMAKE_OSX_DEPLOYMENT_TARGET=$MACOSX_DEPLOYMENT_TARGET"
fi

cmake $cmake_flags ..
make config=release -j$NPROC
make install
cd ../../../..
echo "[BUILD SUCCESS] SUNDIALS"
