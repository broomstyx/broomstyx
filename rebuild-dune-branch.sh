# MKL root directory
export MKL_ROOT=/opt/intel/mkl

cd ..
rm -rf ./*/build-cmake
./dune-common/bin/dunecontrol --opts=./broomstyx/config.opts all
