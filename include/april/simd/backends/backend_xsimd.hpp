#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <xsimd/xsimd.hpp>


#include "april/simd/packed_concept.hpp"

namespace april::simd::internal::xsimd {
    template<typename T, size_t Width = 0> struct Packed;


    template<typename T, size_t Width = 0>
    struct Mask {
        using native_type = Packed<T, Width>::native_type::batch_bool_type;

        Mask() = default;
        Mask(native_type d) : data(d) {}
        Mask(bool val) : data(val) {}

        operator native_type() const { return data; }
        static constexpr size_t size() { return native_type::size; }

        template<typename U>
        requires (size() == Mask<U, Width>::size())
        operator Mask<U, Width>() const {
            if constexpr (sizeof(T) == sizeof(U)) {
                // Zero-cost hardware reinterpret for same-size types
                return { ::xsimd::bitwise_cast<typename Mask<U, Width>::native_type>(data) };
            } else {
                // Generates hardware widening/narrowing instructions
                return { ::xsimd::batch_bool_cast<U>(data) };
            }
        }


        // DATA LOADS
        static Mask load(const bool* ptr) { return { native_type::load_unaligned(ptr) }; }
        static Mask load_aligned(const bool* ptr) { return { native_type::load_aligned(ptr) }; }
        static Mask load_unaligned(const bool* ptr) { return { native_type::load_unaligned(ptr) }; }

        // DATA STORES
        void store(bool* ptr) const { data.store_unaligned(ptr); }
        void store_aligned(bool* ptr) const { data.store_aligned(ptr); }
        void store_unaligned(bool* ptr) const { data.store_unaligned(ptr); }

        // Logical Reductions
        friend bool all(const Mask& m) { return ::xsimd::all(m.data); }
        friend bool any(const Mask& m) { return ::xsimd::any(m.data); }
        friend bool none(const Mask& m) { return !any(m); }

        // Bitwise/Logical Ops
        friend Mask operator~(const Mask& m) { return { ~m.data }; }
        friend Mask operator!(const Mask& m) { return { !m.data }; }
        friend Mask operator^(const Mask& lhs, const Mask& rhs) { return { lhs.data ^ rhs.data }; }
        friend Mask operator&&(const Mask& lhs, const Mask& rhs) { return { lhs.data && rhs.data }; }
        friend Mask operator&(const Mask& lhs, const Mask& rhs)  { return { lhs.data & rhs.data }; }
        friend Mask operator||(const Mask& lhs, const Mask& rhs) { return { lhs.data || rhs.data }; }
        friend Mask operator|(const Mask& lhs, const Mask& rhs)  { return { lhs.data | rhs.data }; }

        // equality
        friend Mask operator==(const Mask& lhs, const Mask& rhs) { return { lhs.data == rhs.data }; }
        friend Mask operator!=(const Mask& lhs, const Mask& rhs) { return { lhs.data != rhs.data }; }

        // rotates
        template<unsigned K = 1>
        void rotate_right() {
            rotate<K>();
        }

        template<unsigned K = 1>
        void rotate_left() {
            rotate<(size() - (K % size())) % size()>();
        }


        // EXPORTS / DEBUGGING
        [[nodiscard]] uint64_t to_bitmask() const {
            return data.mask();
        }

        static Mask from_bitmask(uint64_t bits) {
            std::array<bool, size()> lanes{};

            for (size_t i = 0; i < size(); ++i) {
                lanes[i] = ((bits >> i) & uint64_t{1}) != 0;
            }

            return load_unaligned(lanes.data());
        }

        [[nodiscard]] std::array<bool, size()> to_array() const {
            alignas(alignof(native_type)) std::array<bool, size()> result;
            store_aligned(result.data());
            return result;
        }

        [[nodiscard]] std::string to_string() const {
            auto arr = to_array();
            std::stringstream ss;
            ss << "[";
            for (size_t i = 0; i < size(); ++i) {
                ss << (arr[i] ? "true" : "false");
                if (i < size() - 1) ss << ", ";
            }
            ss << "]";
            return ss.str();
        }

        native_type data;
    private:
        template<unsigned K>
        void rotate() {
            constexpr unsigned Shift = K % size();

            if constexpr (Shift == 0) {} else if constexpr (
                sizeof(typename native_type::register_type) <
                sizeof(typename native_type::batch_type::register_type)
            ) {
                // AVX-512 stores masks as compact predicate bits.
                constexpr uint64_t ValidBits = [] {
                    if constexpr (size() == 64) return ~uint64_t{0};
                    else return (uint64_t{1} << size()) - 1;
                }();

                const uint64_t bits = data.mask();

                // Right lane rotation corresponds to a left rotation of mask bits.
                data = native_type::from_mask(
                    ((bits << Shift) | (bits >> (size() - Shift))) & ValidBits
                );
            } else {
                // SSE/AVX masks occupy full SIMD lanes, where native rotation is cheaper.
                using Batch = native_type::batch_type;
                data = ::xsimd::rotate_right<Shift>(Batch(data)) != Batch(T{0});
            }
        }
    };

    // Mixed-type bitwise operators
    template <typename T, typename U, size_t Width>
    requires (sizeof(T) == sizeof(U) && !std::is_same_v<T, U>)
    Mask<T, Width> operator&(const Mask<T, Width>& lhs, const Mask<U, Width>& rhs) {
        return lhs & static_cast<Mask<T, Width>>(rhs);
    }

    template <typename T, typename U, size_t Width>
    requires (sizeof(T) == sizeof(U) && !std::is_same_v<T, U>)
    Mask<T, Width> operator|(const Mask<T, Width>& lhs, const Mask<U, Width>& rhs) {
        return lhs | static_cast<Mask<T, Width>>(rhs);
    }

    template <typename T, typename U, size_t Width>
    requires (sizeof(T) == sizeof(U) && !std::is_same_v<T, U>)
    Mask<T, Width> operator^(const Mask<T, Width>& lhs, const Mask<U, Width>& rhs) {
        return lhs ^ static_cast<Mask<T, Width>>(rhs);
    }




    template<typename T, size_t Width>
    struct Packed {
        using value_type = std::remove_cv_t<T>;
        using mask_type = Mask<value_type, Width>;

        using native_type = std::conditional_t<
            Width == 0,
            ::xsimd::batch<value_type>,
            ::xsimd::make_sized_batch_t<value_type, Width>
        >;

        static constexpr size_t size() { return native_type::size; }

        Packed() = default;
        Packed(T scalar) : data(scalar) {}
        Packed(native_type d) : data(d) {}

        Packed& operator=(T scalar) {
            data = native_type(scalar);
            return *this;
        }

        // ----------
        // DATA LOADS
        // ----------
        template<typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        static Packed load(const PtrT* ptr) { return load_unaligned(ptr); }

        template<typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        static Packed load_unaligned(const PtrT* ptr) {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                alignas(alignof(native_type)) value_type temp[size()];
                for (size_t i = 0; i < size(); ++i)
                    temp[i] = static_cast<value_type>(ptr[i]);

                return { native_type::load_aligned(temp) };
            } else {
                return { native_type::load_unaligned(
                    reinterpret_cast<const value_type*>(ptr)
                ) };
            }
        }

        template<typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        static Packed load_aligned(const PtrT* ptr) {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                return load_unaligned(ptr);
            } else {
                return { native_type::load_aligned(
                    reinterpret_cast<const value_type*>(ptr)
                ) };
            }
        }


        // ------------
        // DATA GATHERS
        // ------------
        // Gather via offsets
        template<typename PtrT, typename IndexType>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        static Packed gather(const PtrT* base_addr, const IndexType& offsets) {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                alignas(alignof(native_type)) value_type temp[size()];
                for (size_t i = 0; i < size(); ++i) {
                    temp[i] = static_cast<value_type>(base_addr[offsets.data[i]]);
                }
                return { ::xsimd::load_aligned(temp) };
            } else {
                return { native_type::gather(reinterpret_cast<const value_type*>(base_addr), offsets.data) };
            }
        }

        // Gather via array of pointers
        template<typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        static Packed gather(const PtrT* const* pointers) {
            return gather_impl(pointers, std::make_index_sequence<size()>{});
        }


        // -----------
        // DATA STORES
        // -----------
        // Default store delegates to unaligned
        template <typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        void store(PtrT* ptr) const { store_unaligned(ptr); }

        template <typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        void store_unaligned(PtrT* ptr) const {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                alignas(alignof(native_type)) value_type temp[size()];
                ::xsimd::store_aligned(temp, data);
                for (size_t i = 0; i < size(); ++i) {
                    ptr[i] = static_cast<PtrT>(temp[i]);
                }
            } else {
                ::xsimd::store_unaligned(reinterpret_cast<value_type*>(ptr), data);
            }
        }

        template <typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        void store_aligned(PtrT* ptr) const {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                store_unaligned(ptr);
            } else {
                ::xsimd::store_aligned(reinterpret_cast<value_type*>(ptr), data);
            }
        }

        template<typename PtrT, typename IndexType>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        void scatter(PtrT* base_addr, const IndexType& offsets) const {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                alignas(alignof(native_type)) value_type temp[size()];
                ::xsimd::store_aligned(temp, data);
                for (size_t i = 0; i < size(); ++i) {
                    base_addr[offsets.data[i]] = static_cast<PtrT>(temp[i]);
                }
            } else {
                data.scatter(reinterpret_cast<value_type*>(base_addr), offsets.data);
            }
        }

        // PERMUTES AND SHUFFLES
        template<size_t... Indices>
         [[nodiscard]] Packed permute() const {
            return { ::xsimd::swizzle<Indices...>(data) };
        }
        template<unsigned K = 1>
        [[nodiscard]] Packed rotate_left() const {
            // xsimd only has rotate_right, so we compute the complement
            constexpr unsigned Shift = (size() - (K % size())) % size(); // extra modulo so if K == size we get shift = 0
            return { ::xsimd::rotate_right<Shift>(data) };
        }
        template<unsigned K = 1>
        [[nodiscard]] Packed rotate_right() const {
            return { ::xsimd::rotate_right<K>(data) };
        }

        // ARITHMETIC
        friend Packed operator+(const Packed& rhs) { return { +rhs.data };  }
        friend Packed operator-(const Packed& rhs) { return { -rhs.data };  }
        friend Packed operator+(const Packed& lhs, const Packed& rhs) { return { lhs.data + rhs.data }; }
        friend Packed operator-(const Packed& lhs, const Packed& rhs) { return { lhs.data - rhs.data }; }
        friend Packed operator*(const Packed& lhs, const Packed& rhs) { return { lhs.data * rhs.data }; }
        friend Packed operator/(const Packed& lhs, const Packed& rhs) { return { lhs.data / rhs.data }; }

        Packed& operator+=(const Packed& rhs) { data += rhs.data; return *this; }
        Packed& operator-=(const Packed& rhs) { data -= rhs.data; return *this; }
        Packed& operator*=(const Packed& rhs) { data *= rhs.data; return *this; }
        Packed& operator/=(const Packed& rhs) { data /= rhs.data; return *this; }

        // COMPARISONS
        friend mask_type operator==(const Packed& lhs, const Packed& rhs) { return { lhs.data == rhs.data }; }
        friend mask_type operator!=(const Packed& lhs, const Packed& rhs) { return { lhs.data != rhs.data }; }
        friend mask_type operator<(const Packed& lhs, const Packed& rhs)  { return { lhs.data < rhs.data }; }
        friend mask_type operator<=(const Packed& lhs, const Packed& rhs) { return { lhs.data <= rhs.data }; }
        friend mask_type operator>(const Packed& lhs, const Packed& rhs)  { return { lhs.data > rhs.data }; }
        friend mask_type operator>=(const Packed& lhs, const Packed& rhs) { return { lhs.data >= rhs.data }; }

        // MATH FUNCTIONS
        friend Packed sqrt(const Packed& x) { return { ::xsimd::sqrt(x.data) }; }
        friend Packed rsqrt(const Packed& x) { return { ::xsimd::rsqrt(x.data) }; }
        friend Packed abs(const Packed& x) { return { ::xsimd::abs(x.data) }; }

        // ROUNDING
        friend Packed round(const Packed& x) { return { ::xsimd::round(x.data) }; }
        friend Packed floor(const Packed& x) { return { ::xsimd::floor(x.data) }; }
        friend Packed ceil(const Packed& x)  { return { ::xsimd::ceil(x.data) };  }

        // Min/Max/FMA
        friend Packed min(const Packed& a, const Packed& b) { return { ::xsimd::min(a.data, b.data) }; }
        friend Packed max(const Packed& a, const Packed& b) { return { ::xsimd::max(a.data, b.data) }; }
        friend Packed fma(const Packed& a, const Packed& b, const Packed& c) { return { ::xsimd::fma(a.data, b.data, c.data) }; }

        // BITWISE (strictly constrained to integer types)
        friend Packed operator~(const Packed& rhs) requires std::is_integral_v<T> { return { ~rhs.data }; }
        friend Packed operator&(const Packed& lhs, const Packed& rhs) requires std::is_integral_v<T> { return { lhs.data & rhs.data }; }
        friend Packed operator|(const Packed& lhs, const Packed& rhs) requires std::is_integral_v<T> { return { lhs.data | rhs.data }; }
        friend Packed operator^(const Packed& lhs, const Packed& rhs) requires std::is_integral_v<T> { return { lhs.data ^ rhs.data }; }

        Packed& operator&=(const Packed& rhs) requires std::is_integral_v<T> { data &= rhs.data; return *this; }
        Packed& operator|=(const Packed& rhs) requires std::is_integral_v<T> { data |= rhs.data; return *this; }
        Packed& operator^=(const Packed& rhs) requires std::is_integral_v<T> { data ^= rhs.data; return *this; }

        // REDUCTIONS
        [[nodiscard]] T reduce_add() const { return ::xsimd::reduce_add(data); }
        [[nodiscard]] T reduce_min() const { return ::xsimd::reduce_min(data); }
        [[nodiscard]] T reduce_max() const { return ::xsimd::reduce_max(data); }


        // MASKING
        // Performs: result[i] = mask[i] ? true_val[i] : false_val[i]
        friend Packed select(const mask_type& mask, const Packed& true_value, const Packed& false_value) {
            return { ::xsimd::select(mask.data, true_value.data, false_value.data) };
        }

        // DEBUGGING
        [[nodiscard]] std::array<T, size()> to_array() const {
            alignas(alignof(native_type)) std::array<T, size()> result;
            store_aligned(result.data());
            return result;
        }

        [[nodiscard]] std::string to_string() const {
            std::stringstream ss;
            // Create a temporary buffer on the stack
            alignas(64) T buffer[size()];
            store(buffer); // Uses the existing store_unaligned internally

            ss << "[";
            for (size_t i = 0; i < size(); ++i) {
                ss << buffer[i];
                if (i < size() - 1) ss << ", ";
            }
            ss << "]";
            return ss.str();
        }

    private:
        template<typename PtrT, size_t... Is>
        static Packed gather_impl(const PtrT* const * pointers, std::index_sequence<Is...>) {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                alignas(alignof(native_type)) value_type temp[size()];
                // Fold expression to unpack the pointers
                size_t i = 0;
                ((temp[i++] = static_cast<value_type>(*pointers[Is])), ...);
                return { ::xsimd::load_aligned(temp) };
            } else {
                return { ::xsimd::batch<value_type>(*reinterpret_cast<const value_type*>(pointers[Is])...) };
            }
        }

        native_type data;
    };


    static_assert(IsSimdType<Packed<double>>);
    static_assert(IsSimdMask<Mask<double>>);
    static_assert(IsSimdType<Packed<float>>);
    static_assert(IsSimdType<Packed<uint32_t>>);
    static_assert(IsSimdType<Packed<uint64_t>>);
    static_assert(IsSimdMask<Mask<uint64_t>>);
}


