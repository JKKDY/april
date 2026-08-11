#pragma once

#include <type_traits>

#include "april/simd/packed.hpp"
#include "april/simd/packed_concept.hpp"

namespace april::simd {
    template<typename From, typename To>
    concept StorageCompatible =
    std::same_as<std::remove_cv_t<From>, std::remove_cv_t<To>> ||
    (
        std::is_enum_v<std::remove_cv_t<From>> &&
        std::same_as<std::underlying_type_t<std::remove_cv_t<From>>,std::remove_cv_t<To>> &&
        sizeof(From) == sizeof(To) &&
        alignof(From) == alignof(To)
    );

    template<typename Loc>
    concept IsLocation = requires(Loc loc, const Loc cloc, typename Loc::packed_type packed) {
        typename Loc::storage_type;
        typename Loc::value_type;
        typename Loc::packed_type;

        requires IsSimdType<typename Loc::packed_type>;

        { cloc.load() } -> std::same_as<typename Loc::packed_type>;
    };

    template<typename Loc>
    concept IsWritableLocation =
        IsLocation<Loc> &&
        requires(Loc location, const typename Loc::packed_type value) {
            location.store(value);
        };




    template<typename T, bool = std::is_enum_v<std::remove_cv_t<T>>>
    struct location_value {
        using type = std::remove_cv_t<T>;
    };

    template<typename T>
    struct location_value<T, true> {
        using type = std::underlying_type_t<std::remove_cv_t<T>>;
    };

    template<typename T>
    using location_value_t = location_value<T>::type;




    template<typename T>
    struct ContiguousLocation {
        using value_type =  std::remove_cv_t<T>;
        using storage_type = packed_storage_t<value_type>;
        using packed_type = Packed<value_type>;

        T* ptr;

        [[nodiscard]] packed_type load() const noexcept {
            return packed_type::load(ptr);
        }

        void store(const packed_type& value) const noexcept
            requires (!std::is_const_v<T>)
        {
            value.store(ptr);
        }
    };



    template<typename T>
    ContiguousLocation(T*) -> ContiguousLocation<T>;



    inline constexpr std::ptrdiff_t dynamic_stride = -1;

    template<typename T, std::ptrdiff_t ByteStride = dynamic_stride>
    struct StridedLocation;

    // static stride
    template<typename T, std::ptrdiff_t ByteStride>
    requires (ByteStride != dynamic_stride)
    struct StridedLocation<T, ByteStride> {
        using value_type = std::remove_cv_t<T>;
        using storage_type = packed_storage_t<value_type>;
        using packed_type = Packed<value_type>;

        static_assert(ByteStride > 0, "StridedLocation requires a positive byte stride");
        static_assert(
            ByteStride >= static_cast<std::ptrdiff_t>(sizeof(value_type)),
            "Byte stride cannot be smaller than the field type"
        );

        T* ptr;

        constexpr StridedLocation(T* ptr) noexcept
            : ptr(ptr)
        {}

        constexpr StridedLocation(T& value) noexcept
            : ptr(&value)
        {}

        [[nodiscard]] packed_type load() const noexcept {
            return packed_type::template gather_strided<ByteStride>(ptr);
        }

        void store(const packed_type& value) const noexcept
            requires (!std::is_const_v<T>)
        {
            value.template scatter_strided<ByteStride>(ptr);
        }
    };

    // runtime stride
    template<typename T>
    struct StridedLocation<T, dynamic_stride> {
        using value_type = std::remove_cv_t<T>;
        using storage_type = packed_storage_t<value_type>;
        using packed_type = Packed<value_type>;

        T* ptr;
        std::ptrdiff_t byte_stride;

        constexpr StridedLocation(T* ptr, std::ptrdiff_t byte_stride) noexcept
            : ptr(ptr), byte_stride(byte_stride)
        {}

        constexpr StridedLocation(T& value, std::ptrdiff_t byte_stride) noexcept
            : ptr(&value), byte_stride(byte_stride)
        {}

        [[nodiscard]] packed_type load() const noexcept {
            return packed_type::gather_strided(ptr, byte_stride);
        }

        void store(const packed_type& value) const noexcept
            requires (!std::is_const_v<T>)
        {
            value.scatter_strided(ptr, byte_stride);
        }
    };


    template<typename T>
    StridedLocation(T*, std::ptrdiff_t) -> StridedLocation<T, dynamic_stride>;

    template<typename T>
    StridedLocation(T&, std::ptrdiff_t) -> StridedLocation<T, dynamic_stride>;


    template<std::ptrdiff_t ByteStride, typename T>
    [[nodiscard]] constexpr auto make_strided_location(T& value) noexcept {
        return StridedLocation<T, ByteStride>{value};
    }

    template<std::ptrdiff_t ByteStride, typename T>
    [[nodiscard]] constexpr auto make_strided_location(T* ptr) noexcept {
        return StridedLocation<T, ByteStride>{ptr};
    }


    template<typename T>
    struct GatherLocation {
        using value_type = std::remove_cv_t<T>;
        using storage_type = packed_storage_t<value_type>;
        using packed_type = Packed<value_type>;
        using offsets_type = ByteOffsets<packed_type::size()>;

        T* ptr;
        offsets_type offsets;

        constexpr GatherLocation(T* ptr, const offsets_type& offsets) noexcept
            : ptr(ptr), offsets(offsets)
        {}

        constexpr GatherLocation(T& value, const offsets_type& offsets) noexcept
            : ptr(&value), offsets(offsets)
        {}

        [[nodiscard]] packed_type load() const noexcept {
            return packed_type::gather(ptr, offsets);
        }

        void store(const packed_type& value) const noexcept
            requires (!std::is_const_v<T>)
        {
            value.scatter(ptr, offsets);
        }
    };

    template<typename T, size_t N>
    GatherLocation(T*, const ByteOffsets<N>&)
        -> GatherLocation<T>;

    template<typename T, size_t N>
    GatherLocation(T&, const ByteOffsets<N>&)
        -> GatherLocation<T>;


    static_assert(IsWritableLocation<ContiguousLocation<float>>);
    static_assert(IsWritableLocation<ContiguousLocation<double>>);
    static_assert(IsWritableLocation<ContiguousLocation<int>>);
    static_assert(IsWritableLocation<ContiguousLocation<ParticleState>>);

    static_assert(IsWritableLocation<StridedLocation<float>>);
    static_assert(IsWritableLocation<StridedLocation<double>>);
    static_assert(IsWritableLocation<StridedLocation<int>>);
    static_assert(IsWritableLocation<StridedLocation<ParticleState>>);

    static_assert(IsWritableLocation<StridedLocation<float, 1>>);
    static_assert(IsWritableLocation<StridedLocation<double, 1>>);
    static_assert(IsWritableLocation<StridedLocation<int, 1>>);
    static_assert(IsWritableLocation<StridedLocation<ParticleState, 1>>);

    static_assert(IsWritableLocation<GatherLocation<float>>);
    static_assert(IsWritableLocation<GatherLocation<double>>);
    static_assert(IsWritableLocation<GatherLocation<int>>);
    static_assert(IsWritableLocation<GatherLocation<ParticleState>>);
}