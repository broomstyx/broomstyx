# MKL root directory
export MKL_ROOT=/opt/intel/mkl

cd ..
rm -rf ./broomstyx/build-cmake
./dune-common/bin/dunecontrol --opts=./broomstyx/config.opts --only=broomstyx all
