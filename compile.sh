#!/bin/bash


#gfortran -O3 -ftree-vectorize -mcmodel=medium -fbounds-check -fbacktrace   voids_v4.1.f90 -o voids_v4.1.x #serial compilation
gfortran -O3 -ffree-line-length-512 -ftree-vectorize -mcmodel=medium  -fopenmp  voids_v4.1.f90 -o voids_v4.1.x #parallel compilation