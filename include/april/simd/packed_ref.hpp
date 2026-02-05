#pragma once
#include <april/simd/simd_traits.hpp>
#include <april/simd/packed.hpp>



#define AP_SIMD_PROXY_COMPOUND(OP) \
    PackedRef& operator OP(const PackedT& val) { \
        (PackedT::load(ptr) OP val).store(ptr); \
        return *this; \
    } \
    PackedRef& operator OP(std::floating_point auto scalar) { \
        (PackedT::load(ptr) OP PackedT(static_cast<value_type>(scalar))).store(ptr); \
        return *this; \
    } \
    PackedRef& operator OP(const PackedRef& other) { \
        (PackedT::load(ptr) OP PackedT::load(other.ptr)).store(ptr); \
        return *this; \
    }

#define AP_SIMD_PROXY_BINARY(OP) \
    friend PackedT operator OP(const PackedRef& lhs, const PackedRef& rhs) { \
        return PackedT(lhs) OP PackedT(rhs); \
    } \
    friend PackedT operator OP(const PackedRef& lhs, std::floating_point auto rhs) { \
        return PackedT(lhs) OP PackedT(static_cast<value_type>(rhs)); \
    } \
    friend PackedT operator OP(std::floating_point auto lhs, const PackedRef& rhs) { \
        return PackedT(static_cast<value_type>(lhs)) OP PackedT(rhs); \
    } \
    friend PackedT operator OP(const PackedRef& lhs, const PackedT& rhs) { \
        return PackedT(lhs) OP rhs; \
    } \
    friend PackedT operator OP(const PackedT& lhs, const PackedRef& rhs) { \
        return lhs OP PackedT(rhs); \
    }

#define AP_SIMD_PROXY_COMPARE(OP) \
    friend auto operator OP(const PackedRef& lhs, const PackedRef& rhs) { \
        return PackedT(lhs) OP PackedT(rhs); \
    } \
    friend auto operator OP(const PackedRef& lhs, std::floating_point auto rhs) { \
        return PackedT(lhs) OP PackedT(static_cast<value_type>(rhs)); \
    } \
    friend auto operator OP(std::floating_point auto lhs, const PackedRef& rhs) { \
        return PackedT(static_cast<value_type>(lhs)) OP PackedT(rhs); \
    }



//----------------------
// PACKED REF DEFINITION
//----------------------
namespace april::simd {
    template <typename T, IsSimdType PackedT = Packed<T>> // second template arg for type injection in tests
    struct PackedRef {
        using value_type = T;
        using mask_type  = decltype(PackedT() == PackedT());

        value_type* ptr = nullptr;

        PackedRef() = default;
        PackedRef(const PackedRef &) = default;

        explicit PackedRef(value_type* p) : ptr(p) {}

        template <typename U> requires std::convertible_to<U, T>
        PackedRef(const PackedRef<U>& other): ptr(other.ptr) {};

        // CONVERSIONS & ASSIGNMENT

        // TODO right now we are just using the default load and store, but later we can template on a bool to aligned/unaligned access
        // Implicit Load
        operator PackedT() const {
            return PackedT::load(ptr);
        }

        // Store Value
        PackedRef& operator=(const PackedT& val) {
            val.store(ptr);
            return *this;
        }

        // Store Scalar (Broadcast)
        PackedRef& operator=(value_type scalar) {
            PackedT(scalar).store(ptr);
            return *this;
        }

        // Copy from Proxy (Value Copy, with self-assignment check)
        PackedRef& operator=(const PackedRef& other) {
            if (*this != other) {
                PackedT::load(other.ptr).store(ptr);
            }
            return *this;
        }

        // UNARY ARITHMETIC
        friend PackedT operator-(const PackedRef& self) {
            return -PackedT(self);
        }
        friend PackedT operator+(const PackedRef& self) {
            return PackedT(self);
        }

        // BINARY ARITHMETIC
        AP_SIMD_PROXY_BINARY(+)
        AP_SIMD_PROXY_BINARY(-)
        AP_SIMD_PROXY_BINARY(*)
        AP_SIMD_PROXY_BINARY(/)

        // COMPOUND ASSIGNMENT
        AP_SIMD_PROXY_COMPOUND(+=)
        AP_SIMD_PROXY_COMPOUND(-=)
        AP_SIMD_PROXY_COMPOUND(*=)
        AP_SIMD_PROXY_COMPOUND(/=)

        // COMPARISONS
        AP_SIMD_PROXY_COMPARE(==)
        AP_SIMD_PROXY_COMPARE(!=)
        AP_SIMD_PROXY_COMPARE(<)
        AP_SIMD_PROXY_COMPARE(<=)
        AP_SIMD_PROXY_COMPARE(>)
        AP_SIMD_PROXY_COMPARE(>=)

        // MATH FUNCTIONS (ADL Forwarding)
        friend PackedT sqrt(const PackedRef& p) {
            return sqrt(PackedT(p));
        }

        friend PackedT rsqrt(const PackedRef& p) {
            return rsqrt(PackedT(p));
        }

        friend PackedT abs(const PackedRef& p) {
            return abs(PackedT(p));
        }

        friend PackedT min(const PackedRef& a, const PackedRef& b) {
            return min(PackedT(a), PackedT(b));
        }

        friend PackedT max(const PackedRef& a, const PackedRef& b) {
            return max(PackedT(a), PackedT(b));
        }

        friend PackedT fma(const PackedRef& a, const PackedRef& b, const PackedRef& c) {
            return fma(PackedT(a), PackedT(b), PackedT(c));
        }
    };
}


#undef AP_SIMD_PROXY_COMPOUND
#undef AP_SIMD_PROXY_COMPARE
#undef AP_SIMD_PROXY_BINARY
