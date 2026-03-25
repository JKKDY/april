#pragma once


// ----------------
// EXECUTOR DEFAULT
// ----------------
#if !defined(AP_EXECUTOR_BACKEND_OMP) && \
!defined(AP_EXECUTOR_BACKEND_NATIVE_BARRIER) && \
!defined(AP_EXECUTOR_BACKEND_NATIVE_SPIN) && \
!defined(AP_EXECUTOR_BACKEND_SEQUENTIAL)

    #define AP_EXECUTOR_BACKEND_NATIVE_SPIN // default executor backend
#endif



// ------------
// SIMD DEFAULT
// ------------
#if !defined(AP_SIMD_BACKEND_XSIMD) && \
!defined(AP_SIMD_BACKEND_STD_SIMD)

// note: in c++26 with full support of std::simd the default will change from xsimd to std::simd
#define AP_SIMD_BACKEND_XSIMD // default simd backend
#endif

