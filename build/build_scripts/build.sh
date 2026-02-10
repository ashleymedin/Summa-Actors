#!/bin/bash

# If compiling on a Digital Research Alliance of Canada cluster,
# load the following modules:
# module load StdEnv/2023
# module load gcc/12.3
# module load openblas/0.3.24
# module load openmpi/4.1.5
# module load netcdf-fortran/4.6.1

# If compiling on Anvil, load the following modules:
# module load gcc/14.2.0
# module load openmpi/4.1.6
# module load openblas/0.3.17
# module load netcdf-fortran/4.5.3


# If a library cannot be found by Cmake, you can specify the path like so:
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$HOME/Summa-Actors/utils/dependencies/caf/"
export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH:$HOME/Summa-Actors/utils/dependencies/sundials/"


# -----------------------------------
# If compiling V3 use the folowing
# -----------------------------------
# cmake -B ./cmake_build -S ..
# cmake --build ./cmake_build --target all -j

# -----------------------------------
# If compiling V4 without sundials use the following
# -----------------------------------
  
# cmake -B ./cmake_build -S .. -DUSE_V4=ON
# cmake --build ./cmake_build --target all -j

# -----------------------------------
# If compiling V4 with sundials use the following (default)
# -----------------------------------

cmake -B ./cmake_build -S .. -DUSE_SUNDIALS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build ./cmake_build --target all -j
