#pragma once

#include <experimental/simd>
#include <array>
#include <string>

#include "april/simd/simd_traits.hpp"

namespace april::simd::internal::std_simd {

    // Alias for brevity
    namespace stdx = std::experimental;

    template<typename T>
    struct Mask {
        using native_type = stdx::simd_mask<T>;
        native_type data;

        Mask() = default;
        Mask(native_type d) : data(d) {}
        Mask(bool val) : data(val) {}

        operator native_type() const { return data; }

        friend bool all(const Mask& m) {
            return stdx::all_of(m.data);
        }

        friend bool any(const Mask& m) {
            return stdx::any_of(m.data);
        }

        // Logical Not (!)
        friend Mask operator!(const Mask& m) {
            return { !m.data };
        }

        // Note: std::simd uses && and || for element-wise boolean logic
        friend Mask operator&&(const Mask& lhs, const Mask& rhs) { return { lhs.data && rhs.data }; }
        friend Mask operator||(const Mask& lhs, const Mask& rhs) { return { lhs.data || rhs.data }; }

        // Also support bitwise operators just in case
        friend Mask operator&(const Mask& lhs, const Mask& rhs) { return { lhs.data && rhs.data }; }
        friend Mask operator|(const Mask& lhs, const Mask& rhs) { return { lhs.data || rhs.data }; }

        // equality
        friend Mask operator==(const Mask& lhs, const Mask& rhs) { return { lhs.data == rhs.data }; }
        friend Mask operator!=(const Mask& lhs, const Mask& rhs) { return { lhs.data != rhs.data }; }
    };


    // Width == 0: Use Native ABI (Best fit for hardware, e.g. 4 doubles on AVX2)
    // Width > 0:  Use Fixed Size ABI (Compiler manages register spanning, e.g. 16 doubles)
    template<typename T, size_t Width = 0>
    struct Wide {
        using value_type = T;

        using native_type = std::conditional_t<
            Width == 0,
            stdx::simd<T>,                                     // Default/Native ABI
            stdx::simd<T, stdx::simd_abi::fixed_size<Width>>   // Fixed Size ABI
        >;

        static constexpr size_t size() { return native_type::size(); }

        Wide() = default;
        Wide(T scalar) : data(scalar) {}
        Wide(native_type d) : data(d) {}


        // DATA LOADS
        static Wide load(const T* ptr) {
            native_type tmp;
            tmp.copy_from(ptr, stdx::element_aligned);
            return { tmp };
        }

        static Wide load_aligned(const T* ptr) {
            native_type tmp;
            tmp.copy_from(ptr, stdx::vector_aligned);
            return { tmp };
        }

        static Wide load_unaligned(const T* ptr) {
            native_type tmp;
            tmp.copy_from(ptr, stdx::element_aligned);
            return { tmp };
        }
        // Implemented via Generator Constructor:
        // "Construct a SIMD vector where the i-th element is base[offsets[i]]"
        template<typename IndexType>
        static Wide gather(const T* base_addr, const IndexType& offsets) {
            return { native_type([&](size_t i) { return base_addr[offsets.data[i]]; }) };
        }

        static Wide gather(const T* const* pointers) {
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
        [[nodiscard]] Wide permute() const {
            return { native_type([&](size_t i) {
                constexpr std::array<size_t, sizeof...(Indices)> idxs = {Indices...};
                return data[idxs[i]];
            }) };
        }
        template<unsigned K = 1>
        [[nodiscard]] Wide rotate_left() const {
            return { native_type([&](size_t i) {
                return data[(i + K) % size()];
            }) };
        }
        template<unsigned K = 1>
        [[nodiscard]] Wide rotate_right() const {
            return { native_type([&](size_t i) {
               return data[(i + size() - (K % size())) % size()];
            }) };
        }


        // ARITHMETIC
        Wide operator+() {return *this;}
        Wide operator-() {data = -data; return *this;}
        friend Wide operator+(const Wide& lhs, const Wide& rhs) { return { lhs.data + rhs.data }; }
        friend Wide operator-(const Wide& lhs, const Wide& rhs) { return { lhs.data - rhs.data }; }
        friend Wide operator*(const Wide& lhs, const Wide& rhs) { return { lhs.data * rhs.data }; }
        friend Wide operator/(const Wide& lhs, const Wide& rhs) { return { lhs.data / rhs.data }; }

        Wide& operator+=(const Wide& rhs) { data += rhs.data; return *this; }
        Wide& operator-=(const Wide& rhs) { data -= rhs.data; return *this; }
        Wide& operator*=(const Wide& rhs) { data *= rhs.data; return *this; }
        Wide& operator/=(const Wide& rhs) { data /= rhs.data; return *this; }

        // COMPARISONS
        friend Mask<T> operator==(const Wide& lhs, const Wide& rhs) { return { lhs.data == rhs.data }; }
        friend Mask<T> operator!=(const Wide& lhs, const Wide& rhs) { return { lhs.data != rhs.data }; }
        friend Mask<T> operator<(const Wide& lhs, const Wide& rhs)  { return { lhs.data < rhs.data }; }
        friend Mask<T> operator<=(const Wide& lhs, const Wide& rhs) { return { lhs.data <= rhs.data }; }
        friend Mask<T> operator>(const Wide& lhs, const Wide& rhs)  { return { lhs.data > rhs.data }; }
        friend Mask<T> operator>=(const Wide& lhs, const Wide& rhs) { return { lhs.data >= rhs.data }; }


        // MATH FUNCTIONS
        friend Wide sqrt(const Wide& x) {
            return stdx::sqrt(x.data);
        }
        friend Wide rsqrt(const Wide& x) {
            return static_cast<value_type>(1.0) / sqrt(x.data);
        }
        friend Wide abs(const Wide& x) {
            return stdx::abs(x.data);
        }

        // Min/Max/FMA
        friend Wide min(const Wide& a, const Wide& b) { return stdx::min(a.data, b.data) ; }
        friend Wide max(const Wide& a, const Wide& b) { return stdx::max(a.data, b.data) ; }
        friend Wide fma(const Wide& a, const Wide& b, const Wide& c) { return stdx::fma(a.data, b.data, c.data) ; }


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



    template<typename T> Wide<T> sqrt(const Wide<T>& x) {
        using stdx::sqrt;
        return { sqrt(x.data) };
    }

    template<typename T> Wide<T> rsqrt(const Wide<T>& x) {
        // std::simd has no direct rsqrt, fallback to 1.0 / sqrt
        using stdx::sqrt;
        return { Wide<T>(1.0) / sqrt(x.data) };
    }

    template<typename T> Wide<T> abs(const Wide<T>& x) {
        using stdx::abs;
        return { abs(x.data) };
    }

    template<typename T> Wide<T> min(const Wide<T>& a, const Wide<T>& b) {
        using stdx::min;
        return { min(a.data, b.data) };
    }

    template<typename T> Wide<T> max(const Wide<T>& a, const Wide<T>& b) {
        using stdx::max;
        return { max(a.data, b.data) };
    }

    template<typename T> Wide<T> fma(const Wide<T>& a, const Wide<T>& b, const Wide<T>& c) {
        return { std::experimental::fma(a.data, b.data, c.data) };
    }

    static_assert(IsSimdType<Wide<double>>);
    static_assert(IsSimdMask<Mask<double>>);
}
