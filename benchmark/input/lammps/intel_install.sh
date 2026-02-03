cmake ../cmake \
  -D CMAKE_CXX_COMPILER=icpx \
  -D CMAKE_C_COMPILER=icx \
  -D PKG_INTEL=yes \
  -D PKG_OPENMP=yes \
  -D INTEL_ARCH=cpu \
  -D CMAKE_CXX_FLAGS="-O3 -xHost -qopenmp -qopt-zmm-usage=high -ffast-math -fimf-precision=low" \
  -D CMAKE_EXE_LINKER_FLAGS="-ltbbmalloc" \
  -D BUILD_MPI=no
make -j 6