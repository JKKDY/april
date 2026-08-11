#pragma once

// this implementation is for debugging only

#include <array>
#include <cstdint>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cstddef>

#include "april/simd/packed_concept.hpp"

namespace april::simd::internal::scalar {

    template<typename T, size_t Width = 0>
    struct Mask {
        using value_type = std::remove_cv_t<T>;
        static constexpr size_t lane_count = Width == 0 ? size_t{8} : Width;

        using native_type = std::array<bool, lane_count>;
        using mask_type = Mask;

        static_assert(lane_count > 0, "A SIMD mask must contain at least one lane");

        native_type data;

        // Construction
        Mask() : data{} {}

        Mask(bool value) {
            data.fill(value);
        }

        Mask(const native_type& value) : data(value) {}

        Mask(native_type&& value) : data(std::move(value)) {}

        // Native representation access
        operator native_type() const {
            return data;
        }

        // Lane count
        static constexpr size_t size() {
            return lane_count;
        }

        // Conversion between mask element types
        template<typename U>
        requires (size() == Mask<U, Width>::size())
        operator Mask<U, Width>() const {
            return Mask<U, Width>{data};
        }

        // Loads
        static Mask load(const bool* ptr) {
            return load_unaligned(ptr);
        }

        static Mask load_aligned(const bool* ptr) {
            return load_unaligned(ptr);
        }

        static Mask load_unaligned(const bool* ptr) {
            Mask result;
            for (size_t i = 0; i < size(); ++i) result.data[i] = ptr[i];
            return result;
        }

        // Stores
        void store(bool* ptr) const {
            store_unaligned(ptr);
        }

        void store_aligned(bool* ptr) const {
            store_unaligned(ptr);
        }

        void store_unaligned(bool* ptr) const {
            for (size_t i = 0; i < size(); ++i) ptr[i] = data[i];
        }

        // Logical reductions
        [[nodiscard]] static bool all(const Mask& mask) {
            for (const bool lane : mask.data) {
                if (!lane) return false;
            }
            return true;
        }

        [[nodiscard]] static bool any(const Mask& mask) {
            for (const bool lane : mask.data) {
                if (lane) return true;
            }
            return false;
        }

        [[nodiscard]] static bool none(const Mask& mask) {
            return !any(mask);
        }

        // Logical operators
        friend Mask operator!(const Mask& mask) {
            return map(mask, [](const bool value) { return !value; });
        }

        template<typename U>
        friend Mask operator&&(const Mask& lhs, const Mask<U, Width>& rhs) {
            return zip(lhs, rhs, [](const bool a, const bool b) { return a && b; });
        }

        template<typename U>
        friend Mask operator||(const Mask& lhs, const Mask<U, Width>& rhs) {
            return zip(lhs, rhs, [](const bool a, const bool b) { return a || b; });
        }

        // Bitwise operators
        friend Mask operator~(const Mask& mask) {
            return map(mask, [](const bool value) { return !value; });
        }

        template<typename U>
        friend Mask operator&(const Mask& lhs, const Mask<U, Width>& rhs) {
            return zip(lhs, rhs, [](const bool a, const bool b) { return a & b; });
        }

        template<typename U>
        friend Mask operator|(const Mask& lhs, const Mask<U, Width>& rhs) {
            return zip(lhs, rhs, [](const bool a, const bool b) { return a | b; });
        }

        template<typename U>
        friend Mask operator^(const Mask& lhs, const Mask<U, Width>& rhs) {
            return zip(lhs, rhs, [](const bool a, const bool b) { return a != b; });
        }

        // Lane-wise equality
        template<typename U>
        friend Mask operator==(const Mask& lhs, const Mask<U, Width>& rhs) {
            return zip(lhs, rhs, [](const bool a, const bool b) { return a == b; });
        }

        template<typename U>
        friend Mask operator!=(const Mask& lhs, const Mask<U, Width>& rhs) {
            return zip(lhs, rhs, [](const bool a, const bool b) { return a != b; });
        }

        // Compound operators
        template<typename U>
        Mask& operator&=(const Mask<U, Width>& rhs) {
            return *this = *this & rhs;
        }

        template<typename U>
        Mask& operator|=(const Mask<U, Width>& rhs) {
            return *this = *this | rhs;
        }

        template<typename U>
        Mask& operator^=(const Mask<U, Width>& rhs) {
            return *this = *this ^ rhs;
        }


        // Mutating lane rotations
        template<unsigned K = 1>
        void rotate_left() {
            constexpr size_t Shift = K % size();

            if constexpr (Shift != 0) {
                native_type rotated{};
                for (size_t i = 0; i < size(); ++i) rotated[i] = data[(i + Shift) % size()];
                data = rotated;
            }
        }

        template<unsigned K = 1>
        void rotate_right() {
            constexpr size_t Shift = K % size();

            if constexpr (Shift != 0) {
                native_type rotated{};
                for (size_t i = 0; i < size(); ++i) rotated[i] = data[(i + size() - Shift) % size()];
                data = rotated;
            }
        }

        // Bitmask import/export
        [[nodiscard]] uint64_t to_bitmask() const {
            static_assert(size() <= 64, "Mask bit export supports at most 64 lanes");

            uint64_t bits = 0;
            for (size_t i = 0; i < size(); ++i) {
                bits |= static_cast<uint64_t>(data[i]) << i;
            }
            return bits;
        }

        static Mask from_bitmask(const uint64_t bits) {
            static_assert(size() <= 64, "Mask bit import supports at most 64 lanes");

            Mask result;
            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = ((bits >> i) & uint64_t{1}) != 0;
            }
            return result;
        }

        // Debugging and inspection
        [[nodiscard]] std::array<bool, size()> to_array() const {
            return data;
        }

        [[nodiscard]] std::string to_string() const {
            std::stringstream ss;
            ss << "[";

            for (size_t i = 0; i < size(); ++i) {
                ss << (data[i] ? "true" : "false");
                if (i < size() - 1) ss << ", ";
            }

            ss << "]";
            return ss.str();
        }

    private:
        template<typename Fn>
        static Mask map(const Mask& value, Fn&& fn) {
            Mask result;
            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = static_cast<bool>(fn(value.data[i]));
            }
            return result;
        }

        template<typename U, typename Fn>
        static Mask zip(const Mask& lhs, const Mask<U, Width>& rhs, Fn&& fn) {
            Mask result;
            for (size_t i = 0; i < size(); ++i)
                result.data[i] = static_cast<bool>(fn(lhs.data[i], rhs.data[i]));

            return result;
        }
    };



    template<typename T, size_t Width = 0>
    struct Packed {
        using value_type = std::remove_cv_t<T>;
        using storage_type = packed_storage_t<value_type>;

        static constexpr size_t lane_count = Width == 0 ? size_t{8} : Width;

        using native_type = std::array<storage_type, lane_count>;
        using mask_type = Mask<storage_type, Width>;
        using packed_type = Packed;

        static_assert(lane_count > 0, "A SIMD pack must contain at least one lane");

        native_type data;

        // Construction
        Packed() = default;

        Packed(value_type scalar) {
            data.fill(load_scalar(scalar));
        }

        Packed(const native_type& value) : data(value) {}

        Packed(native_type&& value) : data(std::move(value)) {}

        Packed& operator=(value_type scalar) {
            data.fill(load_scalar(scalar));
            return *this;
        }

        // Lane count
        static constexpr size_t size() {
            return lane_count;
        }


        // ----------------
        // Contiguous loads
        // ----------------
        template<typename PtrT>
        requires (
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type))
        )
        static Packed load(const PtrT* ptr) {
            return load_unaligned(ptr);
        }

        template<typename PtrT>
        requires (
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type))
        )
        static Packed load_aligned(const PtrT* ptr) {
            return load_unaligned(ptr);
        }

        template<typename PtrT>
        requires (
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type))
        )
        static Packed load_unaligned(const PtrT* ptr) {
            Packed result;

            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = load_scalar(ptr[i]);
            }

            return result;
        }



        // -----------------
        // Contiguous stores
        // -----------------
        template<typename PtrT>
        requires (
            !std::is_const_v<PtrT> &&
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type))
        )
        void store(PtrT* ptr) const {
            store_unaligned(ptr);
        }

        template<typename PtrT>
        requires (
            !std::is_const_v<PtrT> &&
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type))
        )
        void store_aligned(PtrT* ptr) const {
            store_unaligned(ptr);
        }

        template<typename PtrT>
        requires (
            !std::is_const_v<PtrT> &&
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type))
        )
        void store_unaligned(PtrT* ptr) const {
            for (size_t i = 0; i < size(); ++i) {
                ptr[i] = store_scalar<PtrT>(data[i]);
            }
        }


        // -------
        // Gather
        // -------

        // Compile-time byte stride
        template<std::ptrdiff_t ByteStride, typename PtrT>
        requires (
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type))
        )
        static Packed gather_strided(const PtrT* ptr) {
            if constexpr (ByteStride == sizeof(PtrT)) {
                return load(ptr);
            } else {
                Packed result;
                const auto* base = reinterpret_cast<const std::byte*>(ptr);

                for (size_t i = 0; i < size(); ++i) {
                    const auto* address = reinterpret_cast<const PtrT*>(
                        base + static_cast<std::ptrdiff_t>(i) * ByteStride
                    );

                    result.data[i] = load_scalar(*address);
                }

                return result;
            }
        }

        // Runtime byte stride
        template<typename PtrT>
        requires (
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type))
        )
        static Packed gather_strided(
            const PtrT* ptr,
            const std::ptrdiff_t byte_stride
        ) {
            if (byte_stride == sizeof(PtrT)) {
                return load(ptr);
            }

            Packed result;
            const auto* base = reinterpret_cast<const std::byte*>(ptr);

            for (size_t i = 0; i < size(); ++i) {
                const auto* address = reinterpret_cast<const PtrT*>(
                    base + static_cast<std::ptrdiff_t>(i) * byte_stride
                );

                result.data[i] = load_scalar(*address);
            }

            return result;
        }

        // Arbitrary byte offsets
        template<typename PtrT, size_t N>
        requires (
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type)) &&
            (N == size())
        )
        static Packed gather(
            const PtrT* ptr,
            const ByteOffsets<N>& offsets
        ) {
            Packed result;
            const auto* base = reinterpret_cast<const std::byte*>(ptr);

            for (size_t i = 0; i < size(); ++i) {
                const auto* address = reinterpret_cast<const PtrT*>(
                    base + offsets.values[i]
                );

                result.data[i] = load_scalar(*address);
            }

            return result;
        }

        // Arbitrary pointers
        template<typename PointerContainer>
        requires requires(const PointerContainer& pointers) {
            pointers[size_t{}];
            requires std::is_pointer_v<
                std::remove_cvref_t<decltype(pointers[size_t{}])>
            >;
        }
        static Packed gather(const PointerContainer& pointers) {
            using pointer_type = std::remove_cvref_t<
                decltype(pointers[size_t{}])
            >;

            using PtrT = std::remove_pointer_t<pointer_type>;
            using raw_type = std::remove_cv_t<PtrT>;

            static_assert(
                std::is_arithmetic_v<raw_type> ||
                std::is_enum_v<raw_type>
            );
            static_assert(sizeof(raw_type) <= sizeof(storage_type));

            Packed result;

            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = load_scalar(*pointers[i]);
            }

            return result;
        }


        // -------
        // Scatter
        // -------

        // Compile-time byte stride
        template<std::ptrdiff_t ByteStride, typename PtrT>
        requires (
            !std::is_const_v<PtrT> &&
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type))
        )
        void scatter_strided(PtrT* ptr) const {
            if constexpr (ByteStride == sizeof(PtrT)) {
                store(ptr);
            } else {
                auto* base = reinterpret_cast<std::byte*>(ptr);

                for (size_t i = 0; i < size(); ++i) {
                    auto* address = reinterpret_cast<PtrT*>(
                        base + static_cast<std::ptrdiff_t>(i) * ByteStride
                    );

                    *address = store_scalar<PtrT>(data[i]);
                }
            }
        }

        // Runtime byte stride
        template<typename PtrT>
        requires (
            !std::is_const_v<PtrT> &&
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type))
        )
        void scatter_strided(
            PtrT* ptr,
            const std::ptrdiff_t byte_stride
        ) const {
            if (byte_stride == sizeof(PtrT)) {
                store(ptr);
                return;
            }

            auto* base = reinterpret_cast<std::byte*>(ptr);

            for (size_t i = 0; i < size(); ++i) {
                auto* address = reinterpret_cast<PtrT*>(
                    base + static_cast<std::ptrdiff_t>(i) * byte_stride
                );

                *address = store_scalar<PtrT>(data[i]);
            }
        }

        // Arbitrary byte offsets
        template<typename PtrT, size_t N>
        requires (
            !std::is_const_v<PtrT> &&
            (std::is_arithmetic_v<std::remove_cv_t<PtrT>> || std::is_enum_v<std::remove_cv_t<PtrT>>) &&
            (sizeof(PtrT) <= sizeof(storage_type)) &&
            (N == size())
        )
        void scatter(
            PtrT* ptr,
            const ByteOffsets<N>& offsets
        ) const {
            auto* base = reinterpret_cast<std::byte*>(ptr);

            for (size_t i = 0; i < size(); ++i) {
                auto* address = reinterpret_cast<PtrT*>(
                    base + offsets.values[i]
                );

                *address = store_scalar<PtrT>(data[i]);
            }
        }

        // Arbitrary pointers
        template<typename PointerContainer>
        requires requires(const PointerContainer& pointers) {
            pointers[size_t{}];
            requires std::is_pointer_v<
                std::remove_cvref_t<decltype(pointers[size_t{}])>
            >;
        }
        void scatter(const PointerContainer& pointers) const {
            using pointer_type = std::remove_cvref_t<
                decltype(pointers[size_t{}])
            >;

            using PtrT = std::remove_pointer_t<pointer_type>;
            using raw_type = std::remove_cv_t<PtrT>;

            static_assert(
                std::is_arithmetic_v<raw_type> ||
                std::is_enum_v<raw_type>
            );
            static_assert(!std::is_const_v<PtrT>);
            static_assert(sizeof(raw_type) <= sizeof(storage_type));

            for (size_t i = 0; i < size(); ++i) {
                *pointers[i] = store_scalar<raw_type>(data[i]);
            }
        }



        // ----------------------
        // Permutes and rotations
        // ----------------------

        template<size_t... Indices>
        [[nodiscard]] Packed permute() const {
            static_assert(sizeof...(Indices) > 0, "A permutation requires at least one index");
            static_assert(((Indices < size()) && ...), "Permutation index is outside the pack");

            constexpr std::array<size_t, sizeof...(Indices)> indices{Indices...};
            Packed result;

            if constexpr (sizeof...(Indices) == 1) {
                result.data.fill(data[indices[0]]);
            } else {
                static_assert(
                    sizeof...(Indices) == size(),
                    "A permutation must provide one index or exactly one index per lane"
                );

                for (size_t i = 0; i < size(); ++i) result.data[i] = data[indices[i]];
            }

            return result;
        }

        template<unsigned K = 1>
        [[nodiscard]] Packed rotate_left() const {
            constexpr size_t Shift = K % size();
            Packed result;

            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = data[(i + Shift) % size()];
            }

            return result;
        }

        template<unsigned K = 1>
        [[nodiscard]] Packed rotate_right() const {
            constexpr size_t Shift = K % size();
            Packed result;

            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = data[(i + size() - Shift) % size()];
            }

            return result;
        }

        // --------------------
        // Arithmetic operators
        // --------------------
        friend Packed operator+(const Packed& value) {
            return value.map([](storage_type lane) {
                return +lane;
            });
        }

        friend Packed operator-(const Packed& value) {
            return value.map([](storage_type lane) {
                return -lane;
            });
        }

        friend Packed operator+(const Packed& lhs, const Packed& rhs) {
            return lhs.zip(rhs, [](storage_type a, storage_type b) {
                return a + b;
            });
        }

        friend Packed operator-(const Packed& lhs, const Packed& rhs) {
            return lhs.zip(rhs, [](storage_type a, storage_type b) {
                return a - b;
            });
        }

        friend Packed operator*(const Packed& lhs, const Packed& rhs) {
            return lhs.zip(rhs, [](storage_type a, storage_type b) {
                return a * b;
            });
        }

        friend Packed operator/(const Packed& lhs, const Packed& rhs) {
            return lhs.zip(rhs, [](storage_type a, storage_type b) {
                return a / b;
            });
        }

        friend Packed operator+(const Packed& lhs, storage_type rhs) {
            return lhs.map([rhs](storage_type a) {
                return a + rhs;
            });
        }

        friend Packed operator+(storage_type lhs, const Packed& rhs) {
            return rhs + lhs;
        }

        friend Packed operator-(const Packed& lhs, storage_type rhs) {
            return lhs.map([rhs](storage_type a) {
                return a - rhs;
            });
        }

        friend Packed operator-(storage_type lhs, const Packed& rhs) {
            return rhs.map([lhs](storage_type b) {
                return lhs - b;
            });
        }

        friend Packed operator*(const Packed& lhs, storage_type rhs) {
            return lhs.map([rhs](storage_type a) {
                return a * rhs;
            });
        }

        friend Packed operator*(storage_type lhs, const Packed& rhs) {
            return rhs * lhs;
        }

        friend Packed operator/(const Packed& lhs, storage_type rhs) {
            return lhs.map([rhs](storage_type a) {
                return a / rhs;
            });
        }

        friend Packed operator/(storage_type lhs, const Packed& rhs) {
            return rhs.map([lhs](storage_type b) {
                return lhs / b;
            });
        }

        // ---------
        // COMPOUNDS
        // ---------
        Packed& operator+=(const Packed& rhs) {
            for (size_t i = 0; i < size(); ++i) data[i] += rhs.data[i];
            return *this;
        }

        Packed& operator-=(const Packed& rhs) {
            for (size_t i = 0; i < size(); ++i) data[i] -= rhs.data[i];
            return *this;
        }

        Packed& operator*=(const Packed& rhs) {
            for (size_t i = 0; i < size(); ++i) data[i] *= rhs.data[i];
            return *this;
        }

        Packed& operator/=(const Packed& rhs) {
            for (size_t i = 0; i < size(); ++i) data[i] /= rhs.data[i];
            return *this;
        }

        Packed& operator+=(storage_type rhs) {
            for (auto& lane : data)
                lane += rhs;
            return *this;
        }

        Packed& operator-=(storage_type rhs) {
            for (auto& lane : data)
                lane -= rhs;
            return *this;
        }

        Packed& operator*=(storage_type rhs) {
            for (auto& lane : data)
                lane *= rhs;
            return *this;
        }

        Packed& operator/=(storage_type rhs) {
            for (auto& lane : data)
                lane /= rhs;
            return *this;
        }

        // -----------
        // Comparisons
        // -----------

        friend mask_type operator==(const Packed& lhs, const Packed& rhs) {
            return lhs.compare(rhs, [](storage_type a, storage_type b) {
                return a == b;
            });
        }

        friend mask_type operator!=(const Packed& lhs, const Packed& rhs) {
            return lhs.compare(rhs, [](storage_type a, storage_type b) {
                return a != b;
            });
        }

        friend mask_type operator<(const Packed& lhs, const Packed& rhs) {
            return lhs.compare(rhs, [](storage_type a, storage_type b) {
                return a < b;
            });
        }

        friend mask_type operator<=(const Packed& lhs, const Packed& rhs) {
            return lhs.compare(rhs, [](storage_type a, storage_type b) {
                return a <= b;
            });
        }

        friend mask_type operator>(const Packed& lhs, const Packed& rhs) {
            return lhs.compare(rhs, [](storage_type a, storage_type b) {
                return a > b;
            });
        }

        friend mask_type operator>=(const Packed& lhs, const Packed& rhs) {
            return lhs.compare(rhs, [](storage_type a, storage_type b) {
                return a >= b;
            });
        }

        friend mask_type operator==(const Packed& lhs, storage_type rhs) {
            return lhs.test([rhs](storage_type a) {
                return a == rhs;
            });
        }

        friend mask_type operator==(storage_type lhs, const Packed& rhs) {
            return rhs == lhs;
        }

        friend mask_type operator!=(const Packed& lhs, storage_type rhs) {
            return lhs.test([rhs](storage_type a) {
                return a != rhs;
            });
        }

        friend mask_type operator!=(storage_type lhs, const Packed& rhs) {
            return rhs != lhs;
        }

        friend mask_type operator<(const Packed& lhs, storage_type rhs) {
            return lhs.test([rhs](storage_type a) {
                return a < rhs;
            });
        }

        friend mask_type operator<(storage_type lhs, const Packed& rhs) {
            return rhs.test([lhs](storage_type b) {
                return lhs < b;
            });
        }

        friend mask_type operator<=(const Packed& lhs, storage_type rhs) {
            return lhs.test([rhs](storage_type a) {
                return a <= rhs;
            });
        }

        friend mask_type operator<=(storage_type lhs, const Packed& rhs) {
            return rhs.test([lhs](storage_type b) {
                return lhs <= b;
            });
        }

        friend mask_type operator>(const Packed& lhs, storage_type rhs) {
            return lhs.test([rhs](storage_type a) {
                return a > rhs;
            });
        }

        friend mask_type operator>(storage_type lhs, const Packed& rhs) {
            return rhs.test([lhs](storage_type b) {
                return lhs > b;
            });
        }

        friend mask_type operator>=(const Packed& lhs, storage_type rhs) {
            return lhs.test([rhs](storage_type a) {
                return a >= rhs;
            });
        }

        friend mask_type operator>=(storage_type lhs, const Packed& rhs) {
            return rhs.test([lhs](storage_type b) {
                return lhs >= b;
            });
        }

        // ---------
        // Selection
        // ---------

        [[nodiscard]] static Packed select(
           const mask_type& mask,
           const Packed& true_value,
           const Packed& false_value
        ) {
            Packed result;

            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = mask.data[i]
                    ? true_value.data[i]
                    : false_value.data[i];
            }

            return result;
        }



        // ----------------
        // Basic operations
        // ----------------

        [[nodiscard]] static Packed abs(const Packed& x) {
            if constexpr (std::is_unsigned_v<storage_type>) {
                return x;
            } else {
                return x.map([](storage_type value) {
                    return std::abs(value);
                });
            }
        }

        [[nodiscard]] static Packed min(const Packed& a, const Packed& b) {
            return a.zip(b, [](storage_type lhs, storage_type rhs) {
                return std::min(lhs, rhs);
            });
        }

        [[nodiscard]] static Packed max(const Packed& a, const Packed& b) {
            return a.zip(b, [](storage_type lhs, storage_type rhs) {
                return std::max(lhs, rhs);
            });
        }

        [[nodiscard]] static Packed clamp(const Packed& x, const Packed& lower, const Packed& upper) {
            return x.zip3(lower, upper, [](storage_type value, storage_type lo, storage_type hi) {
                return std::clamp(value, lo, hi);
            });
        }

        // ----------------
        // Roots and powers
        // ----------------

        [[nodiscard]] static Packed sqrt(const Packed& x) requires std::floating_point<storage_type> {
            return x.map([](storage_type value) {
                return std::sqrt(value);
            });
        }

        [[nodiscard]] static Packed rsqrt(const Packed& x) requires std::floating_point<storage_type> {
            return x.map([](storage_type value) {
                return storage_type{1} / std::sqrt(value);
            });
        }

        [[nodiscard]] static Packed cbrt(const Packed& x) requires std::floating_point<storage_type> {
            return x.map([](storage_type value) {
                return std::cbrt(value);
            });
        }

        [[nodiscard]] static Packed hypot(const Packed& x, const Packed& y) requires std::floating_point<storage_type> {
            return x.zip(y, [](storage_type lhs, storage_type rhs) {
                return std::hypot(lhs, rhs);
            });
        }

        [[nodiscard]] static Packed pow(const Packed& x, const Packed& y)
            requires std::floating_point<storage_type>
        {
            return x.zip(y, [](storage_type lhs, storage_type rhs) {
                return std::pow(lhs, rhs);
            });
        }

        // ---------------------------
        // Exponential and logarithmic
        // ---------------------------

        [[nodiscard]] static Packed exp(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::exp(value);
            });
        }

        [[nodiscard]] static Packed exp2(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::exp2(value);
            });
        }

        [[nodiscard]] static Packed expm1(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::expm1(value);
            });
        }

        [[nodiscard]] static Packed log(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::log(value);
            });
        }

        [[nodiscard]] static Packed ln(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return log(x);
        }

        [[nodiscard]] static Packed log2(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::log2(value);
            });
        }

        [[nodiscard]] static Packed log10(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::log10(value);
            });
        }

        [[nodiscard]] static Packed log1p(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::log1p(value);
            });
        }

        // -------------
        // Trigonometric
        // -------------

        [[nodiscard]] static Packed sin(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::sin(value);
            });
        }

        [[nodiscard]] static Packed cos(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::cos(value);
            });
        }

        [[nodiscard]] static std::pair<Packed, Packed> sincos(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return {sin(x), cos(x)};
        }

        [[nodiscard]] static Packed tan(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::tan(value);
            });
        }

        [[nodiscard]] static Packed asin(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::asin(value);
            });
        }

        [[nodiscard]] static Packed acos(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::acos(value);
            });
        }

        [[nodiscard]] static Packed atan(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::atan(value);
            });
        }

        [[nodiscard]] static Packed atan2(const Packed& y, const Packed& x)
            requires std::floating_point<storage_type>
        {
            return y.zip(x, [](storage_type y_value, storage_type x_value) {
                return std::atan2(y_value, x_value);
            });
        }

        // ----------
        // Hyperbolic
        // ----------

        [[nodiscard]] static Packed sinh(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::sinh(value);
            });
        }

        [[nodiscard]] static Packed cosh(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::cosh(value);
            });
        }

        [[nodiscard]] static Packed tanh(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::tanh(value);
            });
        }

        [[nodiscard]] static Packed asinh(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::asinh(value);
            });
        }

        [[nodiscard]] static Packed acosh(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::acosh(value);
            });
        }

        [[nodiscard]] static Packed atanh(const Packed& x)
            requires std::floating_point<storage_type>
        {
            return x.map([](storage_type value) {
                return std::atanh(value);
            });
        }

        // --------
        // Rounding
        // --------

        [[nodiscard]] static Packed floor(const Packed& x) {
            if constexpr (std::integral<storage_type>) {
                return x;
            } else {
                return x.map([](storage_type value) {
                    return std::floor(value);
                });
            }
        }

        [[nodiscard]] static Packed ceil(const Packed& x) {
            if constexpr (std::integral<storage_type>) {
                return x;
            } else {
                return x.map([](storage_type value) {
                    return std::ceil(value);
                });
            }
        }

        [[nodiscard]] static Packed round(const Packed& x) {
            if constexpr (std::integral<storage_type>) {
                return x;
            } else {
                return x.map([](storage_type value) {
                    return std::round(value);
                });
            }
        }

        [[nodiscard]] static Packed trunc(const Packed& x) {
            if constexpr (std::integral<storage_type>) {
                return x;
            } else {
                return x.map([](storage_type value) {
                    return std::trunc(value);
                });
            }
        }

        [[nodiscard]] static Packed nearbyint(const Packed& x) {
            if constexpr (std::integral<storage_type>) {
                return x;
            } else {
                return x.map([](storage_type value) {
                    return std::nearbyint(value);
                });
            }
        }

        // ------------------
        // Numeric operations
        // ------------------

        [[nodiscard]] static Packed fma(const Packed& x, const Packed& y, const Packed& z) {
            if constexpr (std::integral<storage_type>) {
                return x.zip3(y, z, [](storage_type a, storage_type b, storage_type c) {
                    return a * b + c;
                });
            } else {
                return x.zip3(y, z, [](storage_type a, storage_type b, storage_type c) {
                    return std::fma(a, b, c);
                });
            }
        }

        [[nodiscard]] static Packed fmod(const Packed& x, const Packed& y) requires std::floating_point<storage_type> {
            return x.zip(y, [](storage_type lhs, storage_type rhs) {
                return std::fmod(lhs, rhs);
            });
        }

        [[nodiscard]] static Packed remainder(const Packed& x, const Packed& y)
            requires std::floating_point<storage_type>
        {
            return x.zip(y, [](storage_type lhs, storage_type rhs) {
                return std::remainder(lhs, rhs);
            });
        }

        [[nodiscard]] static Packed copysign(const Packed& magnitude, const Packed& sign)
            requires std::floating_point<storage_type>
        {
            return magnitude.zip(sign, [](storage_type magnitude_value, storage_type sign_value) {
                return std::copysign(magnitude_value, sign_value);
            });
        }

        // --------------
        // Classification
        // --------------

        [[nodiscard]] static mask_type isnan(const Packed& x) {
            if constexpr (std::floating_point<storage_type>) {
                return x.test([](storage_type value) {
                    return std::isnan(value);
                });
            } else {
                return mask_type{false};
            }
        }

        [[nodiscard]] static mask_type isinf(const Packed& x) {
            if constexpr (std::floating_point<storage_type>) {
                return x.test([](storage_type value) {
                    return std::isinf(value);
                });
            } else {
                return mask_type{false};
            }
        }

        [[nodiscard]] static mask_type isfinite(const Packed& x) {
            if constexpr (std::floating_point<storage_type>) {
                return x.test([](storage_type value) {
                    return std::isfinite(value);
                });
            } else {
                return mask_type{true};
            }
        }

        [[nodiscard]] static mask_type signbit(const Packed& x) {
            if constexpr (std::is_unsigned_v<storage_type>) {
                return mask_type{false};
            } else if constexpr (std::integral<storage_type>) {
                return x.test([](storage_type value) {
                    return value < storage_type{0};
                });
            } else {
                return x.test([](storage_type value) {
                    return std::signbit(value);
                });
            }
        }

        // -----------------
        // Bitwise operators
        // -----------------

        friend Packed operator~(const Packed& value) requires std::integral<storage_type> {
            return value.map([](storage_type lane) {
                return static_cast<storage_type>(~lane);
            });
        }

        friend Packed operator&(const Packed& lhs, const Packed& rhs) requires std::integral<storage_type> {
            return lhs.zip(rhs, [](storage_type a, storage_type b) {
                return static_cast<storage_type>(a & b);
            });
        }

        friend Packed operator|(const Packed& lhs, const Packed& rhs) requires std::integral<storage_type> {
            return lhs.zip(rhs, [](storage_type a, storage_type b) {
                return static_cast<storage_type>(a | b);
            });
        }

        friend Packed operator^(const Packed& lhs, const Packed& rhs) requires std::integral<storage_type> {
            return lhs.zip(rhs, [](storage_type a, storage_type b) {
                return static_cast<storage_type>(a ^ b);
            });
        }

        Packed& operator&=(const Packed& rhs) requires std::integral<storage_type> {
            for (size_t i = 0; i < size(); ++i) data[i] &= rhs.data[i];
            return *this;
        }

        Packed& operator|=(const Packed& rhs) requires std::integral<storage_type> {
            for (size_t i = 0; i < size(); ++i) data[i] |= rhs.data[i];
            return *this;
        }

        Packed& operator^=(const Packed& rhs) requires std::integral<storage_type> {
            for (size_t i = 0; i < size(); ++i) data[i] ^= rhs.data[i];
            return *this;
        }

        // ----------
        // Reductions
        // ----------

        [[nodiscard]] storage_type reduce_add() const {
            storage_type result{};

            for (size_t i = 0; i < size(); ++i) result += data[i];

            return result;
        }

        [[nodiscard]] storage_type reduce_min() const {
            storage_type result = data[0];

            for (size_t i = 1; i < size(); ++i) {
                if (data[i] < result) result = data[i];
            }

            return result;
        }

        [[nodiscard]] storage_type reduce_max() const {
            storage_type result = data[0];

            for (size_t i = 1; i < size(); ++i) {
                if (data[i] > result) result = data[i];
            }

            return result;
        }

        // ------------------------
        // Debugging and inspection
        // ------------------------

        [[nodiscard]] std::array<value_type, size()> to_array() const {
            std::array<value_type, size()> result;

            for (size_t i = 0; i < size(); ++i) {
                result[i] = store_scalar<value_type>(data[i]);
            }

            return result;
        }

        [[nodiscard]] std::string to_string() const {
            std::stringstream ss;
            ss << "[";

            for (size_t i = 0; i < size(); ++i) {
                ss << data[i];
                if (i < size() - 1) ss << ", ";
            }

            ss << "]";
            return ss.str();
        }

    private:
        template<typename Fn>
[[nodiscard]] Packed map(Fn&& fn) const {
            Packed result;

            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = static_cast<storage_type>(fn(data[i]));
            }

            return result;
        }

        template<typename Fn>
        [[nodiscard]] Packed zip(const Packed& rhs, Fn&& fn) const {
            Packed result;

            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = static_cast<storage_type>(fn(data[i], rhs.data[i]));
            }

            return result;
        }

        template<typename Fn>
        [[nodiscard]] mask_type compare(const Packed& rhs, Fn&& fn) const {
            mask_type result;

            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = static_cast<bool>(fn(data[i], rhs.data[i]));
            }

            return result;
        }

        template<typename Fn>
        [[nodiscard]] Packed zip3(const Packed& b, const Packed& c, Fn&& fn) const {
            Packed result;

            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = static_cast<storage_type>(fn(data[i], b.data[i], c.data[i]));
            }

            return result;
        }

        template<typename Fn>
        [[nodiscard]] mask_type test(Fn&& fn) const {
            mask_type result;

            for (size_t i = 0; i < size(); ++i) {
                result.data[i] = static_cast<bool>(fn(data[i]));
            }

            return result;
        }

        template<typename U>
        [[nodiscard]] static constexpr storage_type load_scalar(const U value) noexcept {
            if constexpr (std::is_enum_v<U>)
                return static_cast<storage_type>(std::to_underlying(value));
            else
                return static_cast<storage_type>(value);
        }

        template<typename U>
        [[nodiscard]] static constexpr U store_scalar(const storage_type value) noexcept {
            return static_cast<U>(value);
        }
    };

    static_assert(april::simd::IsSimdType<Packed<double>>);
    static_assert(april::simd::IsSimdType<Packed<float>>);
    static_assert(april::simd::IsSimdType<Packed<size_t>>);
    static_assert(april::simd::IsSimdType<Packed<int>>);

    static_assert(april::simd::IsSimdMask<Mask<double>>);
    static_assert(april::simd::IsSimdMask<Mask<float>>);
    static_assert(april::simd::IsSimdMask<Mask<size_t>>);
    static_assert(april::simd::IsSimdMask<Mask<int>>);

    static_assert(april::simd::IsSimdType<Packed<float, 4>>);
    static_assert(april::simd::IsSimdType<Packed<double, 2>>);
    static_assert(april::simd::IsSimdMask<Mask<float, 4>>);
    static_assert(april::simd::IsSimdMask<Mask<double, 2>>);


}