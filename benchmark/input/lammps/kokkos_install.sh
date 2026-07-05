# Clean previous build attempt
rm -rf build
mkdir build && cd build

# The "Maximum SIMD" Configuration for Clang-18
cmake ../cmake \
  -D CMAKE_CXX_COMPILER=clang++-18 \
  -D CMAKE_C_COMPILER=clang-18 \
  -D BUILD_MPI=yes \
  -D PKG_KOKKOS=yes \
  -D Kokkos_ENABLE_SERIAL=yes \
  -D Kokkos_ENABLE_OPENMP=yes \
  -D Kokkos_ARCH_TGL=ON \
  -D CMAKE_CXX_FLAGS="-O3 -march=native -ffast-math" \
  -D CMAKE_BUILD_TYPE=Release

make -j 6