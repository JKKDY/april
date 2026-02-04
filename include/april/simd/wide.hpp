#pragma once
#include <april/simd/simd_traits.hpp>




#if defined(AP_SIMD_BACKEND_XSIMD)
#include "april/simd/backend_xsimd.hpp"

namespace april::simd {
    template<typename T, size_t W = 0> using Wide = internal::xsimd::Wide<T, W>;
    template<typename T> using Mask = internal::xsimd::Mask<T>;
}

#elif  defined(AP_SIMD_BACKEND_STD_SIMD)

#include "april/simd/backend_std_simd.hpp"

namespace april::simd {
    template<typename T, size_t W = 0> using Wide = internal::std_simd::Wide<T, W>;
    template<typename T> using Mask = internal::std_simd::Mask<T>;
}

#endif


#define AP_SIMD_PROXY_COMPOUND(OP) \
    SimdRef& operator OP(const WideT& val) { \
        (WideT::load(ptr) OP val).store(ptr); \
        return *this; \
    } \
    SimdRef& operator OP(std::floating_point auto scalar) { \
        (WideT::load(ptr) OP WideT(static_cast<value_type>(scalar))).store(ptr); \
        return *this; \
    } \
    SimdRef& operator OP(const SimdRef& other) { \
        (WideT::load(ptr) OP WideT::load(other.ptr)).store(ptr); \
        return *this; \
    }

#define AP_SIMD_PROXY_BINARY(OP) \
    friend WideT operator OP(const SimdRef& lhs, const SimdRef& rhs) { \
        return WideT(lhs) OP WideT(rhs); \
    } \
    friend WideT operator OP(const SimdRef& lhs, std::floating_point auto rhs) { \
        return WideT(lhs) OP WideT(static_cast<value_type>(rhs)); \
    } \
    friend WideT operator OP(std::floating_point auto lhs, const SimdRef& rhs) { \
        return WideT(static_cast<value_type>(lhs)) OP WideT(rhs); \
    } \
    friend WideT operator OP(const SimdRef& lhs, const WideT& rhs) { \
        return WideT(lhs) OP rhs; \
    } \
    friend WideT operator OP(const WideT& lhs, const SimdRef& rhs) { \
        return lhs OP WideT(rhs); \
    }

#define AP_SIMD_PROXY_COMPARE(OP) \
    friend auto operator OP(const SimdRef& lhs, const SimdRef& rhs) { \
        return WideT(lhs) OP WideT(rhs); \
    } \
    friend auto operator OP(const SimdRef& lhs, std::floating_point auto rhs) { \
        return WideT(lhs) OP WideT(static_cast<value_type>(rhs)); \
    } \
    friend auto operator OP(std::floating_point auto lhs, const SimdRef& rhs) { \
        return WideT(static_cast<value_type>(lhs)) OP WideT(rhs); \
    }



namespace april::simd {

    template <IsSimdType WideT>
    struct SimdRef {
        using value_type = WideT::value_type;
        using mask_type  = decltype(WideT() == WideT());

        value_type* ptr = nullptr;

        SimdRef() = default;
        explicit SimdRef(value_type* p) : ptr(p) {}

        // CONVERSIONS & ASSIGNMENT

        // Implicit Load
        operator WideT() const {
            return WideT::load(ptr);
        }

        // Store Value
        SimdRef& operator=(const WideT& val) {
            val.store(ptr);
            return *this;
        }

        // Store Scalar (Broadcast)
        SimdRef& operator=(value_type scalar) {
            WideT(scalar).store(ptr);
            return *this;
        }

        // Copy from Proxy (Value Copy, with self-assignment check)
        SimdRef& operator=(const SimdRef& other) {
            if (*this != other) {
                WideT::load(other.ptr).store(ptr);
            }
            return *this;
        }

        // UNARY ARITHMETIC
        friend WideT operator-(const SimdRef& self) {
            return -WideT(self); // Load and negate
        }
        friend WideT operator+(const SimdRef& self) {
            return WideT(self); // Load and return
        }

        // BINARY ARITHMETIC (e.g. proxy + proxy or scalar + proxy) returns Wide
        AP_SIMD_PROXY_BINARY(+)
        AP_SIMD_PROXY_BINARY(-)
        AP_SIMD_PROXY_BINARY(*)
        AP_SIMD_PROXY_BINARY(/)


        // COMPOUND ASSIGNMENT (Proxy op = Value/Scalar)
        AP_SIMD_PROXY_COMPOUND(+=)
        AP_SIMD_PROXY_COMPOUND(-=)
        AP_SIMD_PROXY_COMPOUND(*=)
        AP_SIMD_PROXY_COMPOUND(/=)

        // COMPARISONS (Must return Mask, not bool!)
        // "p.mass == p.density" should behave like "Wide == Wide"
        AP_SIMD_PROXY_COMPARE(==)
        AP_SIMD_PROXY_COMPARE(!=)
        AP_SIMD_PROXY_COMPARE(<)
        AP_SIMD_PROXY_COMPARE(<=)
        AP_SIMD_PROXY_COMPARE(>)
        AP_SIMD_PROXY_COMPARE(>=)


        // MATH FUNCTIONS (ADL Forwarding)
        friend WideT sqrt(const SimdRef& p) {
             return sqrt(WideT(p));
        }

        friend WideT rsqrt(const SimdRef& p) {
             return rsqrt(WideT(p));
        }

        friend WideT abs(const SimdRef& p) {
             return abs(WideT(p));
        }

        friend WideT min(const SimdRef& a, const SimdRef& b) {
             return min(WideT(a), WideT(b));
        }

        friend WideT max(const SimdRef& a, const SimdRef& b) {
             return max(WideT(a), WideT(b));
        }

        friend WideT fma(const SimdRef& a, const SimdRef& b, const SimdRef& c) {
             return fma(WideT(a), WideT(b), WideT(c));
        }
    };
}


#undef AP_SIMD_PROXY_COMPOUND
#undef AP_SIMD_PROXY_COMPARE
#undef AP_SIMD_PROXY_BINARY
