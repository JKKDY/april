#pragma once


// ----------------
// EXECUTOR DEFAULT
// ----------------
#if !defined(APRIL_EXECUTOR_BACKEND_OMP) && \
!defined(APRIL_EXECUTOR_BACKEND_NATIVE_BARRIER) && \
!defined(APRIL_EXECUTOR_BACKEND_NATIVE_SPIN) && \
!defined(APRIL_EXECUTOR_BACKEND_SEQUENTIAL)

    #define APRIL_EXECUTOR_BACKEND_NATIVE_SPIN // default executor backend
#endif



// ------------
// SIMD DEFAULT
// ------------
#if !defined(APRIL_SIMD_BACKEND_XSIMD) && \
!defined(APRIL_SIMD_BACKEND_STD_SIMD)

// note: in c++26 with full support of std::simd the default will change from xsimd to std::simd
#define APRIL_SIMD_BACKEND_XSIMD // default simd backend
#endif

