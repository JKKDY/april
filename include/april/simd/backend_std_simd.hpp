#pragma once

#if __has_include(<simd>)
    #include <simd>
    namespace stdx = std;
#elif __has_include(<experimental/simd>)
    #include <experimental/simd>
    namespace stdx = std::experimental;
#else
    #error "No std::simd support found in the standard library"
#endif
#include <array>


#include <string>

#include "april/simd/simd_traits.hpp"

namespace april::simd::internal::std_simd {

    // Alias for brevity
    #if __has_include(<simd>)
        namespace stdx = std;
    #elif __has_include(<experimental/simd>)
        namespace stdx = std::experimental;
    #else
        #error "No std::simd support found in the standard library"
    #endif

    template<typename T>
    struct Mask {
        using native_type = stdx::simd_mask<T>;
        native_type data;

        Mask() = default;
        Mask(native_type d) : data(d) {}
        Mask(bool val) : data(val) {}

        template<typename U>
        requires (sizeof(T) == sizeof(U))
        operator Mask<U>() const {
            typename Mask<U>::native_type converted;
            // Compilers should optimize this fixed-size loop into a zero-cost cast
            for (size_t i = 0; i < size(); ++i) {
                converted[i] = data[i];
            }
            return { converted };
        }

        operator native_type() const { return data; }
        static constexpr size_t size() { return native_type::size(); }

        // ---------------------
        // DATA LOADS
        // ---------------------
        static Mask load(const bool* ptr) {
            return load_unaligned(ptr);
        }
        static Mask load_aligned(const bool* ptr) {
            Mask m;
            m.data.copy_from(ptr, stdx::vector_aligned);
            return m;
        }
        static Mask load_unaligned(const bool* ptr) {
            Mask m;
            m.data.copy_from(ptr, stdx::element_aligned);
            return m;
        }

        // DATA STORES
        void store(bool* ptr) const { store_unaligned(ptr); }
        void store_aligned(bool* ptr) const {data.copy_to(ptr, stdx::vector_aligned); }
        void store_unaligned(bool* ptr) const { data.copy_to(ptr, stdx::element_aligned); }

        // Logical Reductions
        friend bool all(const Mask& m) { return stdx::all_of(m.data); }
        friend bool any(const Mask& m) { return stdx::any_of(m.data); }
        friend bool none(const Mask& m) { return !any(m); }

        // Bitwise/Logical Ops
        friend Mask operator~(const Mask& m) { return { ~m.data }; }
        friend Mask operator!(const Mask& m) { return { !m.data }; }
        friend Mask operator^(const Mask& lhs, const Mask& rhs) { return { lhs.data ^ rhs.data }; }
        friend Mask operator&&(const Mask& lhs, const Mask& rhs) { return { lhs.data && rhs.data }; }
        friend Mask operator&(const Mask& lhs, const Mask& rhs) { return { lhs.data & rhs.data }; }
        friend Mask operator||(const Mask& lhs, const Mask& rhs) { return { lhs.data || rhs.data }; }
        friend Mask operator|(const Mask& lhs, const Mask& rhs) { return { lhs.data | rhs.data }; }

        // equality
        friend Mask operator==(const Mask& lhs, const Mask& rhs) { return { lhs.data == rhs.data }; }
        friend Mask operator!=(const Mask& lhs, const Mask& rhs) { return { lhs.data != rhs.data }; }

        // ---------------------
        // EXPORTS / DEBUGGING
        // ---------------------
        [[nodiscard]] std::array<bool, native_type::size()> to_array() const {
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
    };


    // Width == 0: Use Native ABI (Best fit for hardware, e.g. 4 doubles on AVX2)
    // Width > 0:  Use Fixed Size ABI (Compiler manages register spanning, e.g. 16 doubles)
    template<typename T, size_t Width = 0>
    struct Packed {
        using value_type = std::remove_const_t<T>;
        using native_type = std::conditional_t<
            Width == 0,
            stdx::simd<value_type>,                                     // Default/Native ABI
            stdx::simd<value_type, stdx::simd_abi::fixed_size<Width>>   // Fixed Size ABI
        >;

        static constexpr size_t size() { return native_type::size(); }

        Packed() = default;
        Packed(T scalar) : data(scalar) {}
        Packed(native_type d) : data(d) {}

        Packed& operator=(T scalar) {
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
                tmp.copy_from(reinterpret_cast<const value_type*>(ptr), stdx::element_aligned);
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
                tmp.copy_from(reinterpret_cast<const value_type*>(ptr), stdx::vector_aligned);
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
                data.copy_to(reinterpret_cast<value_type*>(ptr), stdx::element_aligned);
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
                data.copy_to(reinterpret_cast<value_type*>(ptr), stdx::vector_aligned);
            }
        }

        template<typename PtrT, typename IndexType>
        requires std::is_arithmetic_v<PtrT> && (sizeof(PtrT) <= sizeof(value_type))
        void scatter(PtrT* base_addr, const IndexType& offsets) const {
            // std::simd has no direct scatter -> use scalarized loop (compilers can auto-vectorize)
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
        friend Mask<T> operator==(const Packed& lhs, const Packed& rhs) { return { lhs.data == rhs.data }; }
        friend Mask<T> operator!=(const Packed& lhs, const Packed& rhs) { return { lhs.data != rhs.data }; }
        friend Mask<T> operator<(const Packed& lhs, const Packed& rhs)  { return { lhs.data < rhs.data }; }
        friend Mask<T> operator<=(const Packed& lhs, const Packed& rhs) { return { lhs.data <= rhs.data }; }
        friend Mask<T> operator>(const Packed& lhs, const Packed& rhs)  { return { lhs.data > rhs.data }; }
        friend Mask<T> operator>=(const Packed& lhs, const Packed& rhs) { return { lhs.data >= rhs.data }; }


        // MATH FUNCTIONS
        friend Packed sqrt(const Packed& x) { return stdx::sqrt(x.data); }
        friend Packed rsqrt(const Packed& x) { return static_cast<value_type>(1.0) / sqrt(x.data); }
        friend Packed abs(const Packed& x) { return stdx::abs(x.data); }

        // Min/Max/FMA
        friend Packed min(const Packed& a, const Packed& b) { return stdx::min(a.data, b.data) ; }
        friend Packed max(const Packed& a, const Packed& b) { return stdx::max(a.data, b.data) ; }
        friend Packed fma(const Packed& a, const Packed& b, const Packed& c) { return stdx::fma(a.data, b.data, c.data) ; }

        // rounding
        friend Packed round(const Packed& x) { return { stdx::round(x.data) }; }
        friend Packed floor(const Packed& x) { return { stdx::floor(x.data) }; }
        friend Packed ceil(const Packed& x)  { return { stdx::ceil(x.data) };  }

        // BITWISE OPS
        friend Packed operator~(const Packed& rhs) requires std::is_integral_v<T> { return { ~rhs.data }; }
        friend Packed operator&(const Packed& lhs, const Packed& rhs) requires std::is_integral_v<T> { return { lhs.data & rhs.data }; }
        friend Packed operator|(const Packed& lhs, const Packed& rhs) requires std::is_integral_v<T> { return { lhs.data | rhs.data }; }
        friend Packed operator^(const Packed& lhs, const Packed& rhs) requires std::is_integral_v<T> { return { lhs.data ^ rhs.data }; }

        Packed& operator&=(const Packed& rhs) requires std::is_integral_v<T> { data &= rhs.data; return *this; }
        Packed& operator|=(const Packed& rhs) requires std::is_integral_v<T> { data |= rhs.data; return *this; }
        Packed& operator^=(const Packed& rhs) requires std::is_integral_v<T> { data ^= rhs.data; return *this; }

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
        friend Packed select(const Mask<T>& m, const Packed& true_val, const Packed& false_val) {
            native_type result = false_val.data;
            stdx::where(m.data, result) = true_val.data;
            return { result };
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
            store(buffer); // Uses the existing copy_to internally

            ss << "[";
            for (size_t i = 0; i < size(); ++i) {
                ss << buffer[i];
                if (i < size() - 1) ss << ", ";
            }
            ss << "]";
            return ss.str();
        }
    private:
        native_type data;
    };



    template<typename T> Packed<T> sqrt(const Packed<T>& x) {
        using stdx::sqrt;
        return { sqrt(x.data) };
    }

    template<typename T> Packed<T> rsqrt(const Packed<T>& x) {
        // std::simd has no direct rsqrt, fallback to 1.0 / sqrt
        using stdx::sqrt;
        return { Packed<T>(1.0) / sqrt(x.data) };
    }

    template<typename T> Packed<T> abs(const Packed<T>& x) {
        using stdx::abs;
        return { abs(x.data) };
    }

    template<typename T> Packed<T> min(const Packed<T>& a, const Packed<T>& b) {
        using stdx::min;
        return { min(a.data, b.data) };
    }

    template<typename T> Packed<T> max(const Packed<T>& a, const Packed<T>& b) {
        using stdx::max;
        return { max(a.data, b.data) };
    }

    template<typename T> Packed<T> fma(const Packed<T>& a, const Packed<T>& b, const Packed<T>& c) {
        return { std::experimental::fma(a.data, b.data, c.data) };
    }

    static_assert(IsSimdType<Packed<double>>);
    static_assert(IsSimdType<Packed<float>>);
    static_assert(IsSimdMask<Mask<double>>);
    static_assert(IsSimdType<Packed<uint32_t>>);
    static_assert(IsSimdType<Packed<uint64_t>>);
    static_assert(IsSimdMask<Mask<uint64_t>>);

}
