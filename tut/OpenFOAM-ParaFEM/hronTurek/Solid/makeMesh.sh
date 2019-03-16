#!/usr/bin/env bash

# Check PARAFEM_HOME is set
if [ -z ${PARAFEM_HOME} ]; then
 echo "** WARNING: PARAFEM_HOME not set"
 exit
fi

# Check if softlinks have been made
if [ -L "mg2d" ]; then
  echo "Symlinks exist, continue"
else
  echo "Creating softlinks"
  # Link files to PARAFEM_HOME
  ln -s ${PARAFEM_HOME}/bin/mg2d 
  ln -s ${PARAFEM_HOME}/bin/pf2ensi
  ln -s ${PARAFEM_HOME}/bin/pf2ensi.geo.awk
  ln -s ${PARAFEM_HOME}/bin/pf2ensi.var.awk
fi
# Script to make ensi mesh from paraFem data

# Create mesh for ParaFEM
./mg2d xx24

# Create mesh for ParaVIEW
./pf2ensi xx24
mkdir -p ensi
mv *.ensi.* ensi/.
