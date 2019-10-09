#!/bin/bash

# MKL root directory
export MKL_ROOT=/opt/intel/mkl

# ViennaCL root directory
export VIENNACL_ROOT=/data/dontshootthepianist/ViennaCL/ViennaCL-1.7.1

rm -rf build
mkdir build
cd build

cmake -DENABLE_OPENMP=ON ..
