#!/usr/bin/env bash

# Script to make ensi mesh from paraFem data

./mg2d xx23

./pf2ensi xx23

mkdir -p ensi

mv *.ensi.* ensi/.
