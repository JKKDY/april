#pragma once
#include "april/simd/packed.hpp"


namespace april::simd {

    template<typename T>
    struct MaskedPacked {
        using S = std::remove_const_t<T>;
        using PackedT = Packed<S>;
        using MaskT = PackedMask<S>;

        explicit MaskedPacked(const PackedT& value, const MaskT& mask) noexcept
           : data(value), mask(mask) {}

        
        MaskedPacked() = delete;
        MaskedPacked(const PackedT&, MaskT&&) = delete;

        MaskedPacked(const MaskedPacked&) = default;
        MaskedPacked(MaskedPacked&&) = default;

        MaskedPacked& operator=(const MaskedPacked&) = delete;
        MaskedPacked& operator=(MaskedPacked&&) = delete;

        // implicit conversions
        [[nodiscard]] operator PackedT() const noexcept {
            return data;
        }

        [[nodiscard]] const PackedT& value() const noexcept {
            return data;
        }

        // store data
        MaskedPacked& operator=(const PackedT& rhs) noexcept {
            data = select(mask, rhs, data);
            return *this;
        }

        template<typename Scalar>
        requires std::is_arithmetic_v<Scalar> || std::is_enum_v<Scalar>
        MaskedPacked& operator=(Scalar rhs) noexcept {
            return *this = PackedT(static_cast<S>(rhs));
        }

        // Apply compound arithmetic only to active lanes.
        MaskedPacked& operator+=(const PackedT& rhs) noexcept {
            data = select(mask, data + rhs, data);
            return *this;
        }

        MaskedPacked& operator-=(const PackedT& rhs) noexcept {
            data = select(mask, data - rhs, data);
            return *this;
        }

        MaskedPacked& operator*=(const PackedT& rhs) noexcept {
            data = select(mask, data * rhs, data);
            return *this;
        }

        MaskedPacked& operator/=(const PackedT& rhs) noexcept {
            data = select(mask, data / rhs, data);
            return *this;
        }

        template<typename Scalar>
        requires std::is_arithmetic_v<Scalar>
        MaskedPacked& operator+=(Scalar rhs) noexcept {
            return *this += PackedT(static_cast<S>(rhs));
        }

        template<typename Scalar>
        requires std::is_arithmetic_v<Scalar>
        MaskedPacked& operator-=(Scalar rhs) noexcept {
            return *this -= PackedT(static_cast<S>(rhs));
        }

        template<typename Scalar>
        requires std::is_arithmetic_v<Scalar>
        MaskedPacked& operator*=(Scalar rhs) noexcept {
            return *this *= PackedT(static_cast<S>(rhs));
        }

        template<typename Scalar>
        requires std::is_arithmetic_v<Scalar>
        MaskedPacked& operator/=(Scalar rhs) noexcept {
            return *this /= PackedT(static_cast<S>(rhs));
        }


        // Arithmetic intentionally returns an ordinary packed value.
        friend PackedT operator+(const MaskedPacked& lhs, const PackedT& rhs) noexcept {
            return lhs.data + rhs;
        }

        friend PackedT operator+(const PackedT& lhs, const MaskedPacked& rhs) noexcept {
            return lhs + rhs.data;
        }

        friend PackedT operator-(const MaskedPacked& lhs, const PackedT& rhs) noexcept {
            return lhs.data - rhs;
        }

        friend PackedT operator-(const PackedT& lhs, const MaskedPacked& rhs) noexcept {
            return lhs - rhs.data;
        }

        friend PackedT operator*(const MaskedPacked& lhs, const PackedT& rhs) noexcept {
            return lhs.data * rhs;
        }

        friend PackedT operator*(const PackedT& lhs, const MaskedPacked& rhs) noexcept {
            return lhs * rhs.data;
        }

        friend PackedT operator/(const MaskedPacked& lhs, const PackedT& rhs) noexcept {
            return lhs.data / rhs;
        }

        friend PackedT operator/(const PackedT& lhs, const MaskedPacked& rhs) noexcept {
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
        const MaskT & mask;

    };

}