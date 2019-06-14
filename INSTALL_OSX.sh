#!/bin/bash

# Install brew
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# Install build tools
brew install cmake
brew install llvm # clang-omp
#brew install gcc # Fortran

# Install dependencies (each of these has an enormous number)
brew install cgal
brew install liblas # sqlite, geos, proj, geotiff, qt, etc.
brew install pdal # pcl, pyqt, vtk, etc.

# Checkout the repo
git clone https://github.com/rskelly/libgeo
cd libgeo
git checkout tags/olaf
mkdir build
cd build
cmake -DBUILD_APPS:BOOL=TRUE -DCMAKE_BUILD_TYPE=Release -DOPENSSL_ROOT_DIR=/usr/local/Cellar/openssl/1.0.2s ..
make
sudo make install

echo "Done"

