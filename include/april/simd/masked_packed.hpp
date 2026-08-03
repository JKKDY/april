#pragma once
#include "april/simd/packed.hpp"


namespace april::simd {

    template<typename T, size_t Width = 0>
    struct MaskedPacked {
        using S = std::remove_const_t<T>;

        using packed_type = Packed<S, Width>;
        using mask_type = PackedMask<S, Width>;
        using value_type = packed_type::value_type;

        explicit MaskedPacked(const packed_type& value, const mask_type& mask) noexcept
           : data(value), mask(mask) {}


        MaskedPacked() = delete;
        MaskedPacked(const packed_type&, mask_type&&) = delete;

        MaskedPacked(const MaskedPacked&) = default;
        MaskedPacked(MaskedPacked&&) = default;

        MaskedPacked& operator=(const MaskedPacked&) = delete;
        MaskedPacked& operator=(MaskedPacked&&) = delete;

        // implicit conversions
        [[nodiscard]] operator packed_type() const noexcept {
            return data;
        }

        [[nodiscard]] const packed_type& value() const noexcept {
            return data;
        }

        // store data
        MaskedPacked& operator=(const packed_type& rhs) noexcept {
            data = select(mask, rhs, data);
            return *this;
        }

        template<typename Scalar>
        requires std::is_arithmetic_v<Scalar> || std::is_enum_v<Scalar>
        MaskedPacked& operator=(Scalar rhs) noexcept {
            return *this = packed_type(static_cast<S>(rhs));
        }

        // Apply compound arithmetic only to active lanes.
        MaskedPacked& operator+=(const packed_type& rhs) noexcept {
            data = select(mask, data + rhs, data);
            return *this;
        }

        MaskedPacked& operator-=(const packed_type& rhs) noexcept {
            data = select(mask, data - rhs, data);
            return *this;
        }

        MaskedPacked& operator*=(const packed_type& rhs) noexcept {
            data = select(mask, data * rhs, data);
            return *this;
        }

        MaskedPacked& operator/=(const packed_type& rhs) noexcept {
            data = select(mask, data / rhs, data);
            return *this;
        }

        template<typename Scalar>
        requires std::is_arithmetic_v<Scalar>
        MaskedPacked& operator+=(Scalar rhs) noexcept {
            return *this += packed_type(static_cast<S>(rhs));
        }

        template<typename Scalar>
        requires std::is_arithmetic_v<Scalar>
        MaskedPacked& operator-=(Scalar rhs) noexcept {
            return *this -= packed_type(static_cast<S>(rhs));
        }

        template<typename Scalar>
        requires std::is_arithmetic_v<Scalar>
        MaskedPacked& operator*=(Scalar rhs) noexcept {
            return *this *= packed_type(static_cast<S>(rhs));
        }

        template<typename Scalar>
        requires std::is_arithmetic_v<Scalar>
        MaskedPacked& operator/=(Scalar rhs) noexcept {
            return *this /= packed_type(static_cast<S>(rhs));
        }


        // Arithmetic intentionally returns an ordinary packed value.
        friend packed_type operator+(const MaskedPacked& lhs, const packed_type& rhs) noexcept {
            return lhs.data + rhs;
        }

        friend packed_type operator+(const packed_type& lhs, const MaskedPacked& rhs) noexcept {
            return lhs + rhs.data;
        }

        friend packed_type operator-(const MaskedPacked& lhs, const packed_type& rhs) noexcept {
            return lhs.data - rhs;
        }

        friend packed_type operator-(const packed_type& lhs, const MaskedPacked& rhs) noexcept {
            return lhs - rhs.data;
        }

        friend packed_type operator*(const MaskedPacked& lhs, const packed_type& rhs) noexcept {
            return lhs.data * rhs;
        }

        friend packed_type operator*(const packed_type& lhs, const MaskedPacked& rhs) noexcept {
            return lhs * rhs.data;
        }

        friend packed_type operator/(const MaskedPacked& lhs, const packed_type& rhs) noexcept {
            return lhs.data / rhs;
        }

        friend packed_type operator/(const packed_type& lhs, const MaskedPacked& rhs) noexcept {
            return lhs / rhs.data;
        }


        // rotations
        template<unsigned K = 1>
        void rotate_left() {
            data = data.template rotate_left<K>();
        }

        template<unsigned K = 1>
        void rotate_right() {
            data = data.template rotate_right<K>();
        }

    private:
        Packed<S> data;
        const mask_type & mask;

    };

}