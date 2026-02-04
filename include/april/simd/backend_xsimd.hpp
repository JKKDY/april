#pragma once
#include <xsimd/xsimd.hpp>
#include <utility>
#include <string>
#include <array>

#include "april/simd/simd_traits.hpp"

namespace april::simd::internal::xsimd {

    template<typename T>
    struct Mask {
        using native_type = ::xsimd::batch_bool<T>;
        native_type data;

        Mask() = default;
        Mask(native_type d) : data(d) {}
        Mask(bool val) : data(val) {}

        operator native_type() const { return data; }

        friend bool all(const Mask& m) {
            return ::xsimd::all(m.data);
        }

        friend bool any(const Mask& m) {
            return ::xsimd::any(m.data);
        }

        // Logical Not (!)
        friend Mask operator!(const Mask& m) {
            return { !m.data };
        }

        // Bitwise/Logical And (&& / &)
        friend Mask operator&&(const Mask& lhs, const Mask& rhs) { return { lhs.data && rhs.data }; }
        friend Mask operator&(const Mask& lhs, const Mask& rhs)  { return { lhs.data & rhs.data }; }

        // Bitwise/Logical Or (|| / |)
        friend Mask operator||(const Mask& lhs, const Mask& rhs) { return { lhs.data || rhs.data }; }
        friend Mask operator|(const Mask& lhs, const Mask& rhs)  { return { lhs.data | rhs.data }; }

        // equality
        friend Mask operator==(const Mask& lhs, const Mask& rhs) { return { lhs.data == rhs.data }; }
        friend Mask operator!=(const Mask& lhs, const Mask& rhs) { return { lhs.data != rhs.data }; }
    };

    template<typename T, size_t Width = 0>
    struct Wide {
        using value_type = T;
        using native_type = std::conditional_t<
            Width == 0,
            ::xsimd::batch<T>,
            ::xsimd::make_sized_batch_t<T, Width>
        >;

        static constexpr size_t size() { return native_type::size; }

        Wide() = default;
        Wide(T scalar) : data(scalar) {}
        Wide(native_type d) : data(d) {}


        // DATA LOADS
        static Wide load(const T* ptr) {
            return { ::xsimd::load_unaligned(ptr) };
        }

        static Wide load_aligned(const T* ptr) {
            return { ::xsimd::load_aligned(ptr) };
        }

        static Wide load_unaligned(const T* ptr) {
            return { ::xsimd::load_unaligned(ptr) };
        }

        template<typename IndexType>
        static Wide gather(const T* base_addr, const IndexType& offsets) {
            return { native_type::gather(base_addr, offsets.data) };
        }

        static Wide gather(const T* const* pointers) {
            return gather_impl(pointers, std::make_index_sequence<size()>{});
        }

        // DATA STORES
        void store(T* ptr) const {
            ::xsimd::store_unaligned(ptr, data);
        }

        void store_aligned(T* ptr) const {
            ::xsimd::store_aligned(ptr, data);
        }

        void store_unaligned(T* ptr) const {
            ::xsimd::store_unaligned(ptr, data);
        }

        template<typename IndexType>
        void scatter(T* base_addr, const IndexType& offsets) const {
            data.scatter(base_addr, offsets.data);
        }

        // PERMUTES AND SHUFFLES
        template<size_t... Indices>
         [[nodiscard]] Wide permute() const {
            return { ::xsimd::swizzle<Indices...>(data) };
        }
        template<unsigned K = 1>
        [[nodiscard]] Wide rotate_left() const {
            // xsimd only has rotate_right, so we compute the complement
            constexpr unsigned Shift = size() - (K % size());
            return { ::xsimd::rotate_right<Shift>(data) };
        }
        template<unsigned K = 1>
        [[nodiscard]] Wide rotate_right() const {
            return { ::xsimd::rotate_right<K>(data) };
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
        friend Wide sqrt(const Wide& x) { return { ::xsimd::sqrt(x.data) }; }
        friend Wide rsqrt(const Wide& x) { return { ::xsimd::rsqrt(x.data) }; }
        friend Wide abs(const Wide& x) { return { ::xsimd::abs(x.data) }; }

        // Min/Max/FMA
        friend Wide min(const Wide& a, const Wide& b) { return { ::xsimd::min(a.data, b.data) }; }
        friend Wide max(const Wide& a, const Wide& b) { return { ::xsimd::max(a.data, b.data) }; }
        friend Wide fma(const Wide& a, const Wide& b, const Wide& c) { return { ::xsimd::fma(a.data, b.data, c.data) }; }


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
        template<size_t... Is>
        static Wide gather_impl(const T* const* pointers, std::index_sequence<Is...>) {
            return Wide(::xsimd::batch<T>(*pointers[Is]...));
        }

        native_type data;
    };


    static_assert(IsSimdType<Wide<double>>);
    static_assert(IsSimdMask<Mask<double>>);

}
