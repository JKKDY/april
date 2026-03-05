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
           load_unaligned(ptr);
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

        // DATA LOADS
        static Packed load(const T* ptr) {
            native_type tmp;
            tmp.copy_from(ptr, stdx::element_aligned);
            return { tmp };
        }

        // load a type narrower than T and extend it to match target width
        template<typename NarrowT>
        requires std::is_arithmetic_v<NarrowT> && (sizeof(NarrowT) < sizeof(T)) && std::is_convertible_v<NarrowT, T>
        static Packed load_narrow(const NarrowT* ptr) {
            // create a batch type of the narrow data that has the EXACT same lane count as our target
            using narrow_simd_t = stdx::simd<NarrowT, stdx::simd_abi::fixed_size<size()>>;

            // load strictly the required number of bytes (prevents page boundary segfaults)
            narrow_simd_t narrow_batch;
            narrow_batch.copy_from(ptr, stdx::element_aligned);

            // cast/Extend up to the target width (float->double, uint8_t->uint64_t, etc.)
            return { native_type([&](size_t i) {
                return static_cast<value_type>(narrow_batch[i]);
            }) };
        }

        static Packed load_aligned(const T* ptr) {
            native_type tmp;
            tmp.copy_from(ptr, stdx::vector_aligned);
            return { tmp };
        }

        static Packed load_unaligned(const T* ptr) {
            native_type tmp;
            tmp.copy_from(ptr, stdx::element_aligned);
            return { tmp };
        }
        // Implemented via Generator Constructor:
        // "Construct a SIMD vector where the i-th element is base[offsets[i]]"
        template<typename IndexType>
        static Packed gather(const T* base_addr, const IndexType& offsets) {
            return { native_type([&](size_t i) { return base_addr[offsets.data[i]]; }) };
        }

        static Packed gather(const T* const* pointers) {
            return { native_type([&](size_t i) { return *pointers[i]; }) };
        }


        // DATA STORES
        void store(T* ptr) const {
            data.copy_to(ptr, stdx::element_aligned);
        }

        void store_aligned(T* ptr) const {
            data.copy_to(ptr, stdx::vector_aligned);
        }

        void store_unaligned(T* ptr) const {
            data.copy_to(ptr, stdx::element_aligned);
        }
        // std::simd has no direct scatter, so we scalarize the loop.
        // Compilers (GCC/Clang) are very good at auto-vectorizing this pattern.
        template<typename IndexType>
        void scatter(T* base_addr, const IndexType& offsets) const {
            for (size_t i = 0; i < size(); ++i) {
                base_addr[offsets.data[i]] = data[i];
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
        friend Packed sqrt(const Packed& x) {
            return stdx::sqrt(x.data);
        }
        friend Packed rsqrt(const Packed& x) {
            return static_cast<value_type>(1.0) / sqrt(x.data);
        }
        friend Packed abs(const Packed& x) {
            return stdx::abs(x.data);
        }

        // Min/Max/FMA
        friend Packed min(const Packed& a, const Packed& b) { return stdx::min(a.data, b.data) ; }
        friend Packed max(const Packed& a, const Packed& b) { return stdx::max(a.data, b.data) ; }
        friend Packed fma(const Packed& a, const Packed& b, const Packed& c) { return stdx::fma(a.data, b.data, c.data) ; }

        // rounding
        friend Packed round(const Packed& x) { return { stdx::round(x.data) }; }
        friend Packed floor(const Packed& x) { return { stdx::floor(x.data) }; }
        friend Packed ceil(const Packed& x)  { return { stdx::ceil(x.data) };  }

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
