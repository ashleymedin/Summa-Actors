#! /bin/bash
INSTALL_DIR=$PWD/sundials

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

wget https://github.com/LLNL/sundials/archive/refs/tags/v7.4.0.tar.gz
tar -xzf v7.4.0.tar.gz
cd sundials-7.4.0
mkdir build
cd build
cmake ../ -DBUILD_FORTRAN_MODULE_INTERFACE=ON \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR
make -j 4
make install