#!/bin/bash

# Set some permissions
sudo chown -R $(whoami) /usr/local/share/zsh /usr/local/share/zsh/site-functionschmod u+w /usr/local/share/zsh /usr/local/share/szh/site-functions

# Install brew
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# Install build tools
brew install cmake
brew install llvm # clang-omp
brew install gcc # Fortran

# Install dependencies (each of these has an enormous number)
brew install cgal
brew install liblas # sqlite, geos, proj, geotiff, qt, etc.
brew install pdal # pcl, pyqt, vtk, etc.

# Checkout the repo
git clone https://github.com/rskelly/libgeo
cd libgeo
git checkout -b tag/olaf1.0
mkdir build
cd build
cmake -DBUILD_APPS:BOOL=TRUE -DCMAKE_BUILD_TYPE=Release -DOPENSSL_ROOT_DIR=/usr/local/Cellar/openssl/1.0.2s ..
make
sudo make install

echo "Done"

exit(0);




# /usr/local/opt is where includes and libs are stored



icu4c is keg-only, which means it was not symlinked into /usr/local,
because macOS provides libicucore.dylib (but nothing else).

If you need to have icu4c first in your PATH run:
  echo 'export PATH="/usr/local/opt/icu4c/bin:$PATH"' >> ~/.bash_profile
  echo 'export PATH="/usr/local/opt/icu4c/sbin:$PATH"' >> ~/.bash_profile

For compilers to find icu4c you may need to set:
  export LDFLAGS="-L/usr/local/opt/icu4c/lib"
  export CPPFLAGS="-I/usr/local/opt/icu4c/include"

