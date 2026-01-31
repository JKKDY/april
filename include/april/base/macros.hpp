#pragma once

/// Cross-compiler force-inline
#if defined(_MSC_VER)
#   define AP_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#   define AP_FORCE_INLINE inline __attribute__((always_inline))
#else
#   define AP_FORCE_INLINE inline
#endif


/// Empty Base Optimization / No Unique Address
#if defined(_MSC_VER)
	#define AP_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#else
	#define AP_NO_UNIQUE_ADDRESS [[no_unique_address]]
#endif


/// Pointer aliasing hint
#if defined(_MSC_VER)
	#define AP_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
	#define AP_RESTRICT __restrict__
#else
	#define AP_RESTRICT
#endif

// Prefetch hint
#ifdef _MSC_VER
	#include <intrin.h>
	// MSVC doesn't support the 'rw' (Read/Write) distinction directly in the basic prefetch
	#define AP_PREFETCH(addr) _mm_prefetch(reinterpret_cast<const char*>(addr), _MM_HINT_T0)
	#define AP_PREFETCH_NTA(addr) _mm_prefetch(reinterpret_cast<const char*>(addr), _MM_HINT_NTA)
#elif defined(__GNUC__) || defined(__clang__)
	// GCC / Clang
	#define AP_PREFETCH(addr) __builtin_prefetch(addr, 0, 3)
	#define AP_PREFETCH_NTA(addr) __builtin_prefetch(addr, 0, 0)
#else
	#define AP_PREFETCH(addr)
	#define AP_PREFETCH_NTA(addr)
#endif

