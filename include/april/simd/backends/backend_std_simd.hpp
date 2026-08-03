#pragma once

#if __has_include(<experimental/simd>)
    #include <experimental/simd>
#else
    #error "std::experimental::simd is not available"
#endif

#ifndef __cpp_lib_experimental_parallel_simd
    #error "The standard library does not provide Parallelism TS SIMD"
#endif


#include <array>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <type_traits>

#include "april/simd/packed_concept.hpp"


namespace april::simd::internal::std_simd {

    namespace stdx = std::experimental;

    template<typename T, size_t Width = 0> struct Packed;


    template<typename T, size_t Width>
    struct native_simd_selector {
        using type = stdx::simd<
            std::remove_cv_t<T>,
            stdx::simd_abi::fixed_size<Width>
        >;
    };

    template<typename T>
    struct native_simd_selector<T, 0> {
        using type = stdx::simd<std::remove_cv_t<T>>;
    };

    template<typename T, size_t Width>
    using native_simd_t = typename native_simd_selector<T, Width>::type;


    template<typename T, size_t Width = 0>
    struct Mask {
        using value_type = std::remove_cv_t<T>;
        using native_type = typename native_simd_t<value_type, Width>::mask_type;

        native_type data;

        Mask() = default;
        Mask(native_type d) : data(d) {}
        Mask(bool val) : data(val) {}

        operator native_type() const { return data; }
        static constexpr size_t size() { return native_type::size(); }

        template<typename U>
        requires (size() == Mask<U, Width>::size())
        operator Mask<U, Width>() const {
            typename Mask<U, Width>::native_type converted;

            // std::simd has no portable mask cast, so copy the logical lanes.
            // The fixed-size loop is compile-time bounded and optimizes well.
            for (size_t i = 0; i < size(); ++i) {
                converted[i] = data[i];
            }

            return { converted };
        }

        // ---------------------
        // DATA LOADS
        // ---------------------
        static Mask load(const bool* ptr) {
            return load_unaligned(ptr);
        }

        static Mask load_aligned(const bool* ptr) {
            return load_unaligned(ptr);
        }

        static Mask load_unaligned(const bool* ptr) {
            Mask m;
            for (size_t i = 0; i < size(); ++i) {
                m.data[i] = ptr[i];
            }
            return m;
        }

        // DATA STORES
        void store(bool* ptr) const {
            store_unaligned(ptr);
        }

        void store_aligned(bool* ptr) const {
            store_unaligned(ptr);
        }

        void store_unaligned(bool* ptr) const {
            for (size_t i = 0; i < size(); ++i) {
                ptr[i] = data[i];
            }
        }

        // Logical Reductions
        friend bool all(const Mask& m) { return stdx::all_of(m.data); }
        friend bool any(const Mask& m) { return stdx::any_of(m.data); }
        friend bool none(const Mask& m) { return !any(m); }

        // Bitwise/Logical Ops
        friend Mask operator~(const Mask& m) { return { !m.data }; }
        friend Mask operator!(const Mask& m) { return { !m.data }; }
        friend Mask operator^(const Mask& lhs, const Mask& rhs) { return { lhs.data ^ rhs.data }; }
        friend Mask operator&&(const Mask& lhs, const Mask& rhs) { return { lhs.data && rhs.data }; }
        friend Mask operator&(const Mask& lhs, const Mask& rhs) { return { lhs.data & rhs.data }; }
        friend Mask operator||(const Mask& lhs, const Mask& rhs) { return { lhs.data || rhs.data }; }
        friend Mask operator|(const Mask& lhs, const Mask& rhs) { return { lhs.data | rhs.data }; }

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

        // ---------------------
        // EXPORTS / DEBUGGING
        // ---------------------
        [[nodiscard]] uint64_t to_bitmask() const {
            static_assert(size() <= 64, "Mask bit export supports at most 64 lanes");

            uint64_t bits = 0;
            for (size_t i = 0; i < size(); ++i) {
                bits |= static_cast<uint64_t>(data[i]) << i;
            }
            return bits;
        }

        static Mask from_bitmask(uint64_t bits) {
            static_assert(size() <= 64, "Mask bit import supports at most 64 lanes");

            native_type result;
            for (size_t i = 0; i < size(); ++i) {
                result[i] = ((bits >> i) & uint64_t{1}) != 0;
            }
            return { result };
        }

        [[nodiscard]] std::array<bool, size()> to_array() const {
            // alignas ensures the array starts at a vector-aligned memory boundary
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

    private:
        template<unsigned K>
        void rotate() {
            constexpr unsigned Shift = K % size();

            if constexpr (Shift != 0) {
                native_type rotated;

                // std::simd has no portable mask shuffle, so remap the lanes.
                for (size_t i = 0; i < size(); ++i) {
                    rotated[i] = data[(i + size() - Shift) % size()];
                }

                data = rotated;
            }
        }
    };


    // Mixed-type bitwise AND
    template <typename T, typename U, size_t Width>
    requires (
        Mask<T, Width>::size() == Mask<U, Width>::size() &&
        !std::is_same_v<T, U>
    )
    Mask<T, Width> operator&(
        const Mask<T, Width>& lhs,
        const Mask<U, Width>& rhs
    ) {
        return lhs & static_cast<Mask<T, Width>>(rhs);
    }

    // Mixed-type bitwise OR
    template <typename T, typename U, size_t Width>
    requires (
        Mask<T, Width>::size() == Mask<U, Width>::size() &&
        !std::is_same_v<T, U>
    )
    Mask<T, Width> operator|(
        const Mask<T, Width>& lhs,
        const Mask<U, Width>& rhs
    ) {
        return lhs | static_cast<Mask<T, Width>>(rhs);
    }

    // Mixed-type bitwise XOR
    template <typename T, typename U, size_t Width>
    requires (
        Mask<T, Width>::size() == Mask<U, Width>::size() &&
        !std::is_same_v<T, U>
    )
    Mask<T, Width> operator^(
        const Mask<T, Width>& lhs,
        const Mask<U, Width>& rhs
    ) {
        return lhs ^ static_cast<Mask<T, Width>>(rhs);
    }


    // Width == 0: Use Native ABI (Best fit for hardware, e.g. 4 doubles on AVX2)
    // Width > 0:  Use Fixed Size ABI (Compiler manages register spanning, e.g. 16 doubles)
    template<typename T, size_t Width>
    struct Packed {
        using value_type = std::remove_cv_t<T>;
        using native_type = native_simd_t<value_type, Width>;
        using mask_type = Mask<value_type, Width>;

        static constexpr size_t size() { return native_type::size(); }

        Packed() = default;
        Packed(value_type scalar) : data(scalar) {}
        Packed(native_type d) : data(d) {}

        Packed& operator=(value_type scalar) {
            data = scalar;
            return *this;
        }

        //-----------
        // DATA LOADS
        //-----------
        template<typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        static Packed load(const PtrT* ptr) {
            return load_unaligned(ptr);
        }

        template<typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        static Packed load_unaligned(const PtrT* ptr) {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                // Generator constructor: safe upcasting from narrow memory.
                // This ensures we only read exactly size() elements, satisfying ASAN.
                return { native_type([&](size_t i) {
                    return static_cast<value_type>(ptr[i]);
                }) };
            } else {
                native_type tmp;
                tmp.copy_from(
                    reinterpret_cast<const value_type*>(ptr),
                    stdx::element_aligned
                );
                return { tmp };
            }
        }

        template<typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        static Packed load_aligned(const PtrT* ptr) {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                // Alignment usually only applies to the native register width.
                // For narrow types, we fall back to the safe generator load.
                return load_unaligned(ptr);
            } else {
                native_type tmp;
                tmp.copy_from(
                    reinterpret_cast<const value_type*>(ptr),
                    stdx::vector_aligned
                );
                return { tmp };
            }
        }

        // ------------
        // DATA GATHERS
        // ------------
        // Gather via offsets: handles upcasting PtrT -> value_type
        template<typename PtrT, typename IndexType>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        static Packed gather(const PtrT* base_addr, const IndexType& offsets) {
            // Implemented via Generator Constructor:
            // "Construct a SIMD vector where the i-th element is base[offsets[i]]"
            return { native_type([&](size_t i) {
                return static_cast<value_type>(base_addr[offsets.data[i]]);
            }) };
        }

        // Gather via array of pointers: handles upcasting PtrT -> value_type
        template<typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        static Packed gather(const PtrT* const* pointers) {
            return { native_type([&](size_t i) {
                return static_cast<value_type>(*pointers[i]);
            }) };
        }

        // -----------
        // DATA STORES
        // -----------
        // Default store delegates to unaligned
        template <typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        void store(PtrT* ptr) const {
            store_unaligned(ptr);
        }

        template <typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        void store_unaligned(PtrT* ptr) const {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                // Safe truncation loop: explicitly bound to size() elements.
                // This prevents ASAN from flagging 'over-writes' on narrow buffers.
                for (size_t i = 0; i < size(); ++i) {
                    ptr[i] = static_cast<PtrT>(data[i]);
                }
            } else {
                data.copy_to(
                    reinterpret_cast<value_type*>(ptr),
                    stdx::element_aligned
                );
            }
        }

        template <typename PtrT>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        void store_aligned(PtrT* ptr) const {
            if constexpr (sizeof(PtrT) < sizeof(value_type)) {
                // Alignment usually only applies to native vector widths.
                // Fall back to safe loop for narrow types.
                store_unaligned(ptr);
            } else {
                data.copy_to(
                    reinterpret_cast<value_type*>(ptr),
                    stdx::vector_aligned
                );
            }
        }

        template<typename PtrT, typename IndexType>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        void scatter(PtrT* base_addr, const IndexType& offsets) const {
            // std::simd has no direct scatter -> use scalarized loop
            for (size_t i = 0; i < size(); ++i) {
                base_addr[offsets.data[i]] = static_cast<PtrT>(data[i]);
            }
        }

        // PERMUTES AND SHUFFLES
        // Uses generator + constexpr array to map compile-time indices to runtime generator access
        template<size_t... Indices>
        [[nodiscard]] Packed permute() const {
            return { native_type([&](size_t i) {
                constexpr std::array<size_t, sizeof...(Indices)> idxs = {Indices...};
                return data[idxs[i]];
            }) };
        }

        template<unsigned K = 1>
        [[nodiscard]] Packed rotate_left() const {
            return { native_type([&](size_t i) {
                return data[(i + K) % size()];
            }) };
        }

        template<unsigned K = 1>
        [[nodiscard]] Packed rotate_right() const {
            return { native_type([&](size_t i) {
                return data[(i + size() - (K % size())) % size()];
            }) };
        }

        // ARITHMETIC
        friend Packed operator+(const Packed& rhs) { return { +rhs.data }; }
        friend Packed operator-(const Packed& rhs) { return { -rhs.data }; }
        friend Packed operator+(const Packed& lhs, const Packed& rhs) { return { lhs.data + rhs.data }; }
        friend Packed operator-(const Packed& lhs, const Packed& rhs) { return { lhs.data - rhs.data }; }
        friend Packed operator*(const Packed& lhs, const Packed& rhs) { return { lhs.data * rhs.data }; }
        friend Packed operator/(const Packed& lhs, const Packed& rhs) { return { lhs.data / rhs.data }; }

        Packed& operator+=(const Packed& rhs) { data += rhs.data; return *this; }
        Packed& operator-=(const Packed& rhs) { data -= rhs.data; return *this; }
        Packed& operator*=(const Packed& rhs) { data *= rhs.data; return *this; }
        Packed& operator/=(const Packed& rhs) { data /= rhs.data; return *this; }

        // COMPARISONS
        friend mask_type operator==(const Packed& lhs, const Packed& rhs) {
            return compare_impl(lhs, rhs, [](const auto& a, const auto& b) {
                return a == b;
            });
        }

        friend mask_type operator!=(const Packed& lhs, const Packed& rhs) {
            return compare_impl(lhs, rhs, [](const auto& a, const auto& b) {
                return a != b;
            });
        }

        friend mask_type operator<(const Packed& lhs, const Packed& rhs) {
            return compare_impl(lhs, rhs, [](const auto& a, const auto& b) {
                return a < b;
            });
        }

        friend mask_type operator<=(const Packed& lhs, const Packed& rhs) {
            return compare_impl(lhs, rhs, [](const auto& a, const auto& b) {
                return a <= b;
            });
        }

        friend mask_type operator>(const Packed& lhs, const Packed& rhs) {
            return compare_impl(lhs, rhs, [](const auto& a, const auto& b) {
                return a > b;
            });
        }

        friend mask_type operator>=(const Packed& lhs, const Packed& rhs) {
            return compare_impl(lhs, rhs, [](const auto& a, const auto& b) {
                return a >= b;
            });
        }

        // MATH FUNCTIONS
        friend Packed sqrt(const Packed& x) {
            return { stdx::sqrt(x.data) };
        }

        friend Packed rsqrt(const Packed& x) {
            return {
                native_type(value_type{1}) / stdx::sqrt(x.data)
            };
        }

        friend Packed abs(const Packed& x) {
            return { stdx::abs(x.data) };
        }

        // Min/Max/FMA
        friend Packed min(const Packed& a, const Packed& b) {
            return { stdx::min(a.data, b.data) };
        }

        friend Packed max(const Packed& a, const Packed& b) {
            return { stdx::max(a.data, b.data) };
        }

        friend Packed fma(const Packed& a, const Packed& b, const Packed& c) {
            return { stdx::fma(a.data, b.data, c.data) };
        }

        // rounding
        friend Packed round(const Packed& x) { return { stdx::round(x.data) }; }
        friend Packed floor(const Packed& x) { return { stdx::floor(x.data) }; }
        friend Packed ceil(const Packed& x)  { return { stdx::ceil(x.data) }; }

        // BITWISE OPS
        friend Packed operator~(const Packed& rhs)
            requires std::is_integral_v<value_type> {
            return { ~rhs.data };
        }

        friend Packed operator&(const Packed& lhs, const Packed& rhs)
            requires std::is_integral_v<value_type> {
            return { lhs.data & rhs.data };
        }

        friend Packed operator|(const Packed& lhs, const Packed& rhs)
            requires std::is_integral_v<value_type> {
            return { lhs.data | rhs.data };
        }

        friend Packed operator^(const Packed& lhs, const Packed& rhs)
            requires std::is_integral_v<value_type> {
            return { lhs.data ^ rhs.data };
        }

        Packed& operator&=(const Packed& rhs)
            requires std::is_integral_v<value_type> {
            data &= rhs.data;
            return *this;
        }

        Packed& operator|=(const Packed& rhs)
            requires std::is_integral_v<value_type> {
            data |= rhs.data;
            return *this;
        }

        Packed& operator^=(const Packed& rhs)
            requires std::is_integral_v<value_type> {
            data ^= rhs.data;
            return *this;
        }

        // REDUCTION
        [[nodiscard]] value_type reduce_add() const {
            // stdx::reduce defaults to std::plus<>() which is transparent
            return stdx::reduce(data);
        }

        [[nodiscard]] value_type reduce_min() const {
            return stdx::reduce(data, [](const auto& a, const auto& b) {
                return stdx::min(a, b);
            });
        }

        [[nodiscard]] value_type reduce_max() const {
            return stdx::reduce(data, [](const auto& a, const auto& b) {
                return stdx::max(a, b);
            });
        }

        // MASKING
        // Performs: result[i] = mask[i] ? true_val[i] : false_val[i]
        friend Packed select(
            const mask_type& m,
            const Packed& true_val,
            const Packed& false_val
        ) {
            native_type result = false_val.data;
            stdx::where(m.data, result) = true_val.data;
            return { result };
        }

        // DEBUGGING
        [[nodiscard]] std::array<value_type, size()> to_array() const {
            alignas(alignof(native_type)) std::array<value_type, size()> result;
            store_aligned(result.data());
            return result;
        }

        [[nodiscard]] std::string to_string() const {
            std::stringstream ss;

            // Create a temporary buffer on the stack
            alignas(64) value_type buffer[size()];
            store(buffer);

            ss << "[";
            for (size_t i = 0; i < size(); ++i) {
                ss << buffer[i];
                if (i < size() - 1) ss << ", ";
            }
            ss << "]";
            return ss.str();
        }

    private:
        template<typename Compare>
        static mask_type compare_impl(
            const Packed& lhs,
            const Packed& rhs,
            Compare compare
        ) {
            // libstdc++ experimental::simd has a broken fixed-size comparison
            // path for unsigned 64-bit types on LP64 platforms.
            if constexpr (Width != 0 && sizeof(value_type) == 8) {
                typename mask_type::native_type result;

                for (size_t i = 0; i < size(); ++i) {
                    result[i] = compare(lhs.data[i], rhs.data[i]);
                }

                return { result };
            } else {
                return { compare(lhs.data, rhs.data) };
            }
        }

        native_type data;
    };


    static_assert(IsSimdType<Packed<double>>);
    static_assert(IsSimdType<Packed<float>>);
    static_assert(IsSimdType<Packed<size_t>>);
    static_assert(IsSimdType<Packed<int>>);

    static_assert(IsSimdMask<Mask<double>>);
    static_assert(IsSimdMask<Mask<float>>);
    static_assert(IsSimdMask<Mask<size_t>>);
    static_assert(IsSimdMask<Mask<int>>);
}
