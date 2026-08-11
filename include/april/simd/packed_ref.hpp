#pragma once
#include "april/simd/packed_concept.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/locations.hpp"

#define AP_SIMD_PROXY_COMPOUND(OP) \
	PackedRef& operator OP(const packed_type& val) { \
		*this = static_cast<packed_type>(*this) OP val; \
		return *this; \
	} \
	template<typename Scalar> requires std::is_arithmetic_v<Scalar> \
	PackedRef& operator OP(Scalar scalar) { \
		*this = static_cast<packed_type>(*this) OP packed_type(static_cast<value_type>(scalar)); \
		return *this; \
	} \
	PackedRef& operator OP(const PackedRef& other) { \
		*this = static_cast<packed_type>(*this) OP static_cast<packed_type>(other); \
		return *this; \
	}

#define AP_SIMD_PROXY_BINARY(OP) \
    friend packed_type operator OP(const PackedRef& lhs, const PackedRef& rhs) { \
        return static_cast<packed_type>(lhs) OP static_cast<packed_type>(rhs); \
    } \
    template <typename Scalar> requires std::is_arithmetic_v<Scalar> \
    friend packed_type operator OP(const PackedRef& lhs, Scalar rhs) { \
        return static_cast<packed_type>(lhs) OP packed_type(static_cast<value_type>(rhs)); \
    } \
    template <typename Scalar> requires std::is_arithmetic_v<Scalar> \
    friend packed_type operator OP(Scalar lhs, const PackedRef& rhs) { \
        return packed_type(static_cast<value_type>(lhs)) OP static_cast<packed_type>(rhs); \
    } \
    friend packed_type operator OP(const PackedRef& lhs, const packed_type& rhs) { \
        return static_cast<packed_type>(lhs) OP rhs; \
    } \
    friend packed_type operator OP(const packed_type& lhs, const PackedRef& rhs) { \
        return lhs OP static_cast<packed_type>(rhs); \
    }

#define AP_SIMD_PROXY_COMPARE(OP) \
	friend auto operator OP(const PackedRef& lhs, const PackedRef& rhs) { \
		return static_cast<packed_type>(lhs) OP static_cast<packed_type>(rhs); \
	} \
	template<april::simd::IsPackableValue Scalar> \
	friend auto operator OP(const PackedRef& lhs, Scalar rhs) { \
		return static_cast<packed_type>(lhs) OP packed_type(static_cast<value_type>(rhs)); \
	} \
	template<april::simd::IsPackableValue Scalar> \
	friend auto operator OP(Scalar lhs, const PackedRef& rhs) { \
		return packed_type(static_cast<value_type>(lhs)) OP static_cast<packed_type>(rhs); \
	}


//----------------------
// PACKED REF DEFINITION
//----------------------
namespace april::simd {


    template<typename T, bool = std::is_enum_v<std::remove_const_t<T>>>
    struct packed_memory {
        using type = std::remove_const_t<T>;
    };

    template<typename T>
    struct packed_memory<T, true> {
        using type = std::underlying_type_t<std::remove_const_t<T>>;
    };

    template<typename T>
    using packed_memory_t = packed_memory<T>::type;

    // T dictates the physical pointer width in memory.
    // packed_type dictates the hardware register width (defaults to matching T).
    // Overriding packed_type allows for loading different data types
    // i.e. trivially convertable or narrower types e.g. loading floats from memory into a simd double register
    template<IsLocation Location>
    struct PackedRef {
        using storage_type = Location::storage_type;
        using value_type = Location::value_type;
        using packed_type = Location::packed_type;
        using mask_type = packed_type::mask_type;

        Location location;

        PackedRef(const PackedRef&) = default;

        explicit PackedRef(Location location)
            : location(std::move(location))
        {}

        template<IsLocation ULocation>
        requires
            std::constructible_from<Location, ULocation> &&
            std::same_as<typename ULocation::packed_type, packed_type>
        PackedRef(const PackedRef<ULocation>& other)
            : location(other.location)
        {}



        // CONVERSIONS & ASSIGNMENT
        // TODO right now we are just using the default load and store, but later we can template on a bool to aligned/unaligned access
        // Implicit Load
        [[nodiscard]] packed_type load() const {
            return location.load();
        }

        operator packed_type() const {
            return load();
        }


        // Packed value
        PackedRef& operator=(const packed_type& value)
            requires IsWritableLocation<Location>
        {
            location.store(value);
            return *this;
        }

        // Semantic scalar
        template<typename Scalar>
        requires (
            IsWritableLocation<Location> &&
            IsPackableValue<std::remove_cvref_t<Scalar>> &&
            std::convertible_to<Scalar, value_type>
        )
        PackedRef& operator=(Scalar scalar)
        {
            location.store(
                packed_type(static_cast<value_type>(scalar))
            );
            return *this;
        }

        // Packed-like/proxy value, e.g. MaskedPacked
        template<typename P>
        requires (
            IsWritableLocation<Location> &&
            !IsPackableValue<std::remove_cvref_t<P>> &&
            std::convertible_to<const P&, packed_type>
        )
        PackedRef& operator=(const P& value)
        {
            location.store(static_cast<packed_type>(value));
            return *this;
        }


        // UNARY ARITHMETIC
        friend packed_type operator-(const PackedRef& self) {
            return -packed_type(self);
        }
        friend packed_type operator+(const PackedRef& self) {
            return packed_type(self);
        }

        // PERMUTES AND SHUFFLES
        template<size_t... Indices>
         [[nodiscard]] packed_type permute() const {
            packed_type val = *this;
            return val.template permute<Indices...>();
        }
        template<unsigned K = 1>
        [[nodiscard]] packed_type rotate_left() const {
            packed_type val = *this;
            return val.template rotate_left<K>();
        }

        template<unsigned K = 1>
        [[nodiscard]] packed_type rotate_right() const {
            packed_type val = *this;
            return val.template rotate_right<K>();
        }

        [[nodiscard]] auto to_array() const {
            return location.load().to_array();
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
    };
}


#undef AP_SIMD_PROXY_COMPOUND
#undef AP_SIMD_PROXY_COMPARE
#undef AP_SIMD_PROXY_BINARY
















