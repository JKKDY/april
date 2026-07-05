#pragma once
#include <april/simd/simd_traits.hpp>
#include <april/simd/packed.hpp>



#define AP_SIMD_PROXY_COMPOUND(OP) \
    PackedRef& operator OP(const PackedT& val) { \
        (static_cast<PackedT>(*this) OP val).store(ptr); \
        return *this; \
    } \
    template <typename Scalar> requires std::is_arithmetic_v<Scalar> \
    PackedRef& operator OP(Scalar scalar) { \
        (static_cast<PackedT>(*this) OP PackedT(static_cast<value_type>(scalar))).store(ptr); \
        return *this; \
    } \
    PackedRef& operator OP(const PackedRef& other) { \
        (static_cast<PackedT>(*this) OP static_cast<PackedT>(other)).store(ptr); \
        return *this; \
    }

#define AP_SIMD_PROXY_BINARY(OP) \
    friend PackedT operator OP(const PackedRef& lhs, const PackedRef& rhs) { \
        return static_cast<PackedT>(lhs) OP static_cast<PackedT>(rhs); \
    } \
    template <typename Scalar> requires std::is_arithmetic_v<Scalar> \
    friend PackedT operator OP(const PackedRef& lhs, Scalar rhs) { \
        return static_cast<PackedT>(lhs) OP PackedT(static_cast<value_type>(rhs)); \
    } \
    template <typename Scalar> requires std::is_arithmetic_v<Scalar> \
    friend PackedT operator OP(Scalar lhs, const PackedRef& rhs) { \
        return PackedT(static_cast<value_type>(lhs)) OP static_cast<PackedT>(rhs); \
    } \
    friend PackedT operator OP(const PackedRef& lhs, const PackedT& rhs) { \
        return static_cast<PackedT>(lhs) OP rhs; \
    } \
    friend PackedT operator OP(const PackedT& lhs, const PackedRef& rhs) { \
        return lhs OP static_cast<PackedT>(rhs); \
    }

#define AP_SIMD_PROXY_COMPARE(OP) \
    friend auto operator OP(const PackedRef& lhs, const PackedRef& rhs) { \
        return static_cast<PackedT>(lhs) OP static_cast<PackedT>(rhs); \
    } \
    template <typename Scalar> requires std::is_arithmetic_v<Scalar> \
    friend auto operator OP(const PackedRef& lhs, Scalar rhs) { \
        return static_cast<PackedT>(lhs) OP PackedT(static_cast<value_type>(rhs)); \
    } \
    template <typename Scalar> requires std::is_arithmetic_v<Scalar> \
    friend auto operator OP(Scalar lhs, const PackedRef& rhs) { \
        return PackedT(static_cast<value_type>(lhs)) OP static_cast<PackedT>(rhs); \
    }


//----------------------
// PACKED REF DEFINITION
//----------------------
namespace april::simd {

    // T dictates the physical pointer width in memory.
    // PackedT dictates the hardware register width (defaults to matching T).
    // Overriding PackedT allows for loading different data types
    // i.e. trivially convertable or narrower types e.g. loading floats from memory into a simd double register
    template <typename T, IsSimdType PackedT = Packed<std::remove_const_t<T>>>
    struct PackedRef {
        using memory_type = std::remove_const_t<T>;
        using value_type = PackedT::value_type;
        using mask_type  = decltype(PackedT() == PackedT());

        using ptr_type = std::conditional_t<std::is_const_v<PackedT>, const T*, T*>;

        ptr_type ptr = nullptr;

        PackedRef() = default;
        PackedRef(const PackedRef &) = default;

        explicit PackedRef(ptr_type p) : ptr(p) {}

        template <typename U_Mem, IsSimdType U_Packed>
        requires std::convertible_to<U_Mem, memory_type>
        PackedRef(const PackedRef<U_Mem, U_Packed>& other): ptr(other.ptr) {}


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

        // Store Scalar (broadcast, handles floats, ints, and strong enums natively)
        template <typename Scalar>
        requires std::is_arithmetic_v<Scalar> || std::is_enum_v<Scalar>
        PackedRef& operator=(Scalar scalar) {
            PackedT(static_cast<value_type>(scalar)).store(ptr);
            return *this;
        }

        // Copy from Proxy (Value Copy, with self-assignment check)
        PackedRef& operator=(const PackedRef& other) {
            if (this->ptr != other.ptr) {
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

        // PERMUTES AND SHUFFLES
        template<size_t... Indices>
         [[nodiscard]] PackedT permute() const {
            PackedT val = *this;
            return val.template permute<Indices...>();
        }
        template<unsigned K = 1>
        [[nodiscard]] PackedT rotate_left() const {
            PackedT val = *this;
            return val.template rotate_left<K>();
        }

        template<unsigned K = 1>
        [[nodiscard]] PackedT rotate_right() const {
            PackedT val = *this;
            return val.template rotate_right<K>();
        }

        [[nodiscard]] auto to_array() const {
            return PackedT::load(ptr).to_array();
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















