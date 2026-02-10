#!/bin/bash
INSTALL_DIR=$PWD/caf

# If compiling on a Digital Research Alliance of Canada cluster,
# load the following modules:
# module load StdEnv/2023
# module load gcc/12.3
# module load openblas/0.3.24
# module load netcdf-fortran/4.6.1

# If compiling on Anvil, load the following modules:
# module load gcc/14.2.0
# module load openblas/0.3.17
# module load netcdf-fortran/4.5.3

wget https://github.com/actor-framework/actor-framework/archive/refs/tags/1.1.0.tar.gz
tar -xvf 1.1.0.tar.gz


echo "Installing CAF to $INSTALL_DIR"

cd actor-framework-1.1.0
./configure --prefix=$INSTALL_DIR
cd build
make -j 8
make install
