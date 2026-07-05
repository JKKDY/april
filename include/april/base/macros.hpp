#pragma once

/// Cross-compiler force-inline
#if defined(_MSC_VER)
#   define APRIL_FORCE_INLINE [[msvc::forceinline]]
#elif defined(__GNUC__) || defined(__clang__)
#   define APRIL_FORCE_INLINE __attribute__((always_inline))
#else
#   define APRIL_FORCE_INLINE
#endif


/// Empty Base Optimization / No Unique Address
#if defined(_MSC_VER)
	#define APRIL_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#else
	#define APRIL_NO_UNIQUE_ADDRESS [[no_unique_address]]
#endif


/// Pointer aliasing hint
#if defined(_MSC_VER)
	#define APRIL_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
	#define APRIL_RESTRICT __restrict__
#else
	#define APRIL_RESTRICT
#endif

// Prefetch hint
#ifdef _MSC_VER
	#include <intrin.h>
	// MSVC doesn't support the 'rw' (Read/Write) distinction directly in the basic prefetch
	#define APRIL_PREFETCH(addr) _mm_prefetch(reinterpret_cast<const char*>(addr), _MM_HINT_T0)
	#define APRIL_PREFETCH_NTA(addr) _mm_prefetch(reinterpret_cast<const char*>(addr), _MM_HINT_NTA)
#elif defined(__GNUC__) || defined(__clang__)
	// GCC / Clang
	#define APRIL_PREFETCH(addr) __builtin_prefetch(addr, 0, 3)
	#define APRIL_PREFETCH_NTA(addr) __builtin_prefetch(addr, 0, 0)
#else
	#define APRIL_PREFETCH(addr)
	#define APRIL_PREFETCH_NTA(addr)
#endif


#if defined(__clang__)
	#define APRIL_PRAGMA(x) _Pragma(#x)
	#define APRIL_UNROLL_LOOP() APRIL_PRAGMA(unroll)
	#define APRIL_UNROLL_LOOP_N(N) APRIL_PRAGMA(unroll N)
#else
	// not clang: no-op
	#define APRIL_UNROLL_LOOP()
	#define APRIL_UNROLL_LOOP_N(N)
#endif













