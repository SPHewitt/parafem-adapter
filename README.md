# parafem-adapter
This repository contains a parafem-adapter for preCICE. The adapter is primarily developed for the coupling of OpenFOAM and paraFEM for fluid-structure interaction problems however the framework that is similar to that of OpenFPCI can be used to develop more advanced adapters.

## ---- Description----
The repository contains a tutorial based on the Turek and Hron benchmark. The tutorial uses OpenFOAM for the fluid solver and ParaFEM for the solid solver.

## ---- Installation ----
The follwing links contain detailed instructions of how to compile [preCICE](https://github.com/precice/precice/wiki) and the [openfoam-adapter](https://github.com/precice/openfoam-adapter) on both linux machines and HPC services. The following sections provide instructions to install [ParaFEM](https://github.com/leemargetts/ParaFEM) and the [parafem-adapter](https://github.com/SPHewitt/parafem-adapter) on Ubuntu 18.04 and full install instructions for the UK's National Computing service ARCHER.

### Ubuntu 18.04
** TO COMPLETE

### ARCHER
To build preCICE, the openfoam-adpater and parafem-adapter on ARCHER, use the following instructions.

Load the following modules:
> module swap PrgEnv-cray PrgEnv-gnu  
> module load cray-petsc   
> module load cmake/3.10.2  
> module load eigen/3.3.0  

Use the following environment variable throught the insall process
> export CRAYPE_LINK_TYPE=dynamic  

#### boost
The boost version on ARCHER is too old to be used, so download and build the boost libraries from source at [boost_1.65.1](https://www.boost.org/users/history/version_1_65_1.html).

1. Download and extract the source code to a directory of your choice on your work partition.
> mkdir dependencies  
> cd dependencies   
> wget https://dl.bintray.com/boostorg/release/1.65.1/source/boost_1_65_1.tar.bz2  
> tar -xvjf boost_1_65_1.tar.bz2
2. Create an install directory for the boost include and lib directories.
> mkdir boost  
3. Enter into source code of boost, compile and install.
> cd boost_1_65_1  
> CC=cc CXX=CC ./booststrap.sh --with-libraries=log,thread,system,filesystem,program_options,test --prefix=/path/to/users/boost/   
> ./b2 install
4. Add the boost environment variables to your bashrc.
> \#boost  
> export BOOST_ROOT=/path/to/users/boost/  
> export LIBRARY_PATH=${BOOST_ROOT}/lib:$LIBRARY_PATH  
> export LD_LIBRARY_PATH=${BOOST_ROOT}/lib:$LD_LIBRARY_PATH  
> export CPLUS_INCLUDE_PATH=${BOOST_ROOT}/include:$CPLUS_INCLUDE_PATH

#### preCICE
preCICE can now be compiled.

1. Download and extract the latest version of preCICE. It has currently be tested with v1.5.0
> wget https://github.com/precice/precice/archive/v1.5.0.tar.gz  
> tar -zxf v1.5.0.tar.gz  
> cd precice-1.5.0  
2. Build and compile the source code. Note, currently not working with python or petsc.
> mkdir build  
> cd build  
> CC=cc CXX=CC cmake -DPYTHON=OFF -DPETSC=OFF -DBUILD_SHARED_LIBS=ON ..  
> make -j 4
3. Add the following environment variables to your .bashrc file.
> \# preCICE  
> export PRECICE_ROOT=path/to/precice-1.5.0  
> export LD_LIBRARY_PATH=${PRECICE_ROOT}/build:$LD_LIBRARY_PATH

#### yaml-cpp
Before the openfoam-adapter can be compiled we need to install yaml-cpp.

1. Download and extract source file to a directory of your choice
> cd dependencies  
> wget https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.5.3.tar.gz  
> tar -zxf yaml-cpp-0.5.3.tar.gz  
> cd yaml-cpp
2. Compile and install.
> mkdir build  
> cd build  
>  CC=cc CXX=CC cmake -DYAML_CPP_BUILD_TOOLS=OFF -DBUILD_SHARED_LIBS=ON ..  
> make -j 4  

#### openfoam-adapter  
1. The openfoam adapter is currently different for each OpenFOAM version so ensure you fetch the OpenFOAM6 branch once you have cloned the repository.
> git clone https://github.com/precice/openfoam-adapter.git  
> cd openfoam-adapter  
> git checkout OpenFOAM6  
2. Load the openfoam/v6 module and compile the adapter
> module load openfoam/v6  
> source /work/y07/y07/cse/openfoam/v6/build64/OpenFOAM-6/etc/bashrc  
> ./Allwmake  

#### parafem-adapter
TODO:
* Create a new build.inc file
* Change makefile to fit this
*
