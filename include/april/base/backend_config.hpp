#pragma once

#ifndef APRIL_FAST_MATH_ENABLED

    #if defined(__FAST_MATH__) || \
    (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0) || \
    (defined(_M_FP_FAST) && _M_FP_FAST)

        #define APRIL_FAST_MATH_ENABLED 1
    #else
        #define APRIL_FAST_MATH_ENABLED 0
    #endif

#endif


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
!defined(APRIL_SIMD_BACKEND_STD_SIMD) && \
!defined(APRIL_SIMD_BACKEND_SCALAR)

    // will be changed with widespread support of std::simd in the future to APRIL_SIMD_BACKEND_STD_SIMD
    #define APRIL_SIMD_BACKEND_XSIMD
#endif

