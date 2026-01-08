// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Julian Deller-Yee
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

