#pragma once

#include <type_traits>

#include "april/particle/attributes.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/packed_concept.hpp"


namespace april {
    enum class MaskPolicy;
}

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

    template<typename T, typename Mask, size_t Width = double_width>
        requires particle::IsParticleAttributes<std::remove_const_t<T>>
    struct PointerPack {
        using value_type = T;

    private:
        std::array<T*, Width> ptrs{};
        size_t offset = 0;
        Mask* mask = nullptr;

    public:
        constexpr PointerPack() = default;

        template<typename U> requires std::convertible_to<U*, T*>
        constexpr explicit PointerPack(U* ptr, Mask& mask) noexcept
        : mask(&mask) {
            for (size_t i = 0; i < Width; ++i)
                ptrs[i] = ptr + i;
        }

        constexpr explicit PointerPack(std::array<T*, Width> ptrs, Mask& mask) noexcept
            : ptrs(ptrs), mask(&mask)
        {}

        template<typename U> requires std::convertible_to<U*, T*>
        constexpr explicit PointerPack(const std::array<U*, Width>& source, Mask& mask) noexcept
        : mask(&mask) {
            for (size_t i = 0; i < Width; ++i)
                ptrs[i] = source[i];
        }

        [[nodiscard]] bool active(const size_t lane) const noexcept {
            const auto bits = mask->to_bitmask();
            return (bits & (decltype(bits){1} << lane)) != 0;
        }

        [[nodiscard]] auto bitmask() const noexcept {
            return mask->to_bitmask();
        }

        [[nodiscard]] auto array_mask() const noexcept {
            return mask->to_array();
        }

        [[nodiscard]] T& operator[](const size_t lane) noexcept {
            return *ptrs[(lane + offset) % Width];
        }

        [[nodiscard]] const T& operator[](size_t lane) const noexcept {
            return *ptrs[(lane + offset) % Width];
        }

        template<unsigned K = 1>
        void rotate_left() noexcept {
            offset = (offset + K) % Width;
        }

        template<unsigned K = 1>
        void rotate_right() noexcept {
            offset = (offset + Width - (K % Width)) % Width;
        }
    };



    enum class Alignment {
        Unaligned,
        Aligned
    };


    template<typename T, Alignment A = Alignment::Unaligned>
    struct ContiguousLocation {
            using value_type = std::remove_cv_t<T>;
            using storage_type = packed_storage_t<value_type>;
            using packed_type = Packed<value_type>;

            T* ptr;

            [[nodiscard]] packed_type load() const noexcept {
                auto* storage = reinterpret_cast<const storage_type*>(ptr);

                if constexpr (A == Alignment::Aligned)
                    return packed_type::load_aligned(storage);
                else
                    return packed_type::load_unaligned(storage);
            }

            void store(const packed_type& value) const noexcept
                requires (!std::is_const_v<T>)
            {
                auto* storage = reinterpret_cast<storage_type*>(ptr);

                if constexpr (A == Alignment::Aligned)
                    value.store_aligned(storage);
                else
                    value.store_unaligned(storage);
            }
        };

    template<typename T>
    ContiguousLocation(T*) -> ContiguousLocation<T>;

    template<typename T>
    auto aligned_location(T* ptr) {
        return ContiguousLocation<T, Alignment::Aligned>{ptr};
    }

    template<typename T>
    auto contiguous_location(T* ptr) {
        return ContiguousLocation<T, Alignment::Unaligned>{ptr};
    }

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

    static_assert(IsWritableLocation<StridedLocation<float, 2 * sizeof(float)>>);
    static_assert(IsWritableLocation<StridedLocation<double, 2 * sizeof(double)>>);
    static_assert(IsWritableLocation<StridedLocation<int, 2 * sizeof(int)>>);
    static_assert(IsWritableLocation<StridedLocation<ParticleState, 2 * sizeof(ParticleState)>>);

    static_assert(IsWritableLocation<GatherLocation<float>>);
    static_assert(IsWritableLocation<GatherLocation<double>>);
    static_assert(IsWritableLocation<GatherLocation<int>>);
    static_assert(IsWritableLocation<GatherLocation<ParticleState>>);
}
