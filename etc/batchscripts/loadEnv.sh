#!/usr/bin/env bash


module swap PrgEnv-cray PrgEnv-gnu
module load cray-petsc
module load cmake/3.10.2
module load eigen/3.3.0

# Build Boost 1.65.1
# TODO: Instructions 

# Build preCICE
# CC=cc CXX=CC cmake -DPYTHON=OFF -DPETSC=OFF -DBUILD_SHARED_LIBS=ON ..
# make -j 4

# Build yaml-cpp 0.5.3
# CC=cc CXX=CC cmake -DYAML_CPP_BUILD_TOOLS=OFF -DBUILD_SHARED_LIBS=ON ..
# make -j 4

# Build OpenFOAM Adapter
# git clone https://github.com/precice/openfoam-adapter.git (OpenFOAM6)
# module load openfoam/v6
# source /work/y07/y07/cse/openfoam/v6/build64/OpenFOAM-6/etc/bashrc
# ./Allwmake 

# Build ParaFEM Adapter
# Build parafem first use gnu compilers
# make 


# HOPE THE TEST WORKS !!!
