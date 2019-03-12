# parafem-adapter
This repository contains a parafem-adapter for preCICE. The adapter is primarily developed for the coupling of OpenFOAM and paraFEM for fluid-structure interaction problems however the framework that is similar to that of OpenFPCI can be used to develop more advanced adapters.

## Description
The repository contains a tutorial based on the Turek and Hron benchmark. The tutorial uses OpenFOAM for the fluid solver and ParaFEM for the solid solver.

## Running
./main precice-config.xml SolverOne MeshOne 
./main precice-config.xml SolverTwo MeshTwo
