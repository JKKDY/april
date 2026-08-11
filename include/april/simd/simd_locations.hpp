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



    template<typename T>
    struct StridedLocation {
        using storage_type = packed_storage_t<T>;
        using value_type = T;
        using packed_type = Packed<storage_type>;

        T* ptr;
        std::ptrdiff_t byte_stride;

        packed_type load() const noexcept;
        void store(const packed_type&) const noexcept
            requires (!std::is_const_v<T>);
    };

    template<typename T>
    StridedLocation(T*, std::ptrdiff_t) -> StridedLocation<T>;



    // template<typename T, IsSimdType OffsetPack>
    // struct GatherLocation {
    //     using storage_type = T;
    //     using value_type = location_value_t<T>;
    //     using packed_type = Packed<value_type>;
    //
    //     T* ptr;
    //     OffsetPack byte_offsets;
    //
    //     packed_type load() const noexcept;
    //     void store(const packed_type&) const noexcept
    //         requires (!std::is_const_v<T>);
    // };
    //
    // template<typename T, size_t Width>
    // GatherLocation(std::array<T*, Width>) -> GatherLocation<T, Width>;



    static_assert(IsLocation<ContiguousLocation<float>>);
    static_assert(IsLocation<ContiguousLocation<double>>);
    static_assert(IsLocation<ContiguousLocation<int>>);

    static_assert(IsLocation<StridedLocation<float>>);
    static_assert(IsLocation<StridedLocation<double>>);
    static_assert(IsLocation<StridedLocation<int>>);

    // static_assert(IsLocation<GatherLocation<float>>);
    // static_assert(IsLocation<GatherLocation<double>>);
    // static_assert(IsLocation<GatherLocation<int>>);

}