#pragma once
#include "april/simd/packed.hpp"
#include "april/utility/debug.hpp"
#include "april/simd/packed_concept.hpp"


namespace april::simd {

    template<
        IsSimdType PackedT,
        IsSimdMask SharedMaskT = typename PackedT::mask_type
    > requires (PackedT::size() == SharedMaskT::size())
    struct MaskedPacked {
        using packed_type = PackedT;
        using mask_type = SharedMaskT;
        using value_type = packed_type::value_type;
        using storage_type = packed_type::storage_type;

        MaskedPacked() = default;
        MaskedPacked(const MaskedPacked&) = default;
        MaskedPacked(MaskedPacked&&) = default;

        MaskedPacked& operator=(const MaskedPacked& rhs) noexcept {
            return *this = rhs.data;
        }

        MaskedPacked& operator=(MaskedPacked&& rhs) noexcept {
            return *this = rhs.data;
        }

        explicit MaskedPacked(const packed_type& value) noexcept
            : data(value), mask(nullptr) {}

        explicit MaskedPacked(const packed_type& value, const mask_type & mask) noexcept
           : data(value), mask(&mask) {}

        void bind_mask(mask_type&&) = delete;
        void bind_mask(mask_type& new_mask) noexcept {
            mask = &new_mask;
        }

        // implicit conversion
        [[nodiscard]] operator packed_type() const noexcept {
            return data;
        }

        [[nodiscard]] const packed_type& value() const noexcept {
            return data;
        }

        // store data
        MaskedPacked& operator=(const packed_type& rhs) noexcept {
            data = select(active_mask(), rhs, data);
            return *this;
        }

        MaskedPacked& operator=(value_type rhs) noexcept {
            return *this = packed_type(rhs);
        }

        // Arithmetic intentionally returns an ordinary packed value.
        [[nodiscard]] packed_type operator+() const noexcept {
            return +data;
        }

        [[nodiscard]] packed_type operator-() const noexcept {
            return -data;
        }

        // Apply compound arithmetic only to active lanes.
        MaskedPacked& operator+=(const packed_type& rhs) noexcept {
            data = select(active_mask(), data + rhs, data);
            return *this;
        }

        MaskedPacked& operator-=(const packed_type& rhs) noexcept {
            data = select(active_mask(), data - rhs, data);
            return *this;
        }

        MaskedPacked& operator*=(const packed_type& rhs) noexcept {
            data = select(active_mask(), data * rhs, data);
            return *this;
        }

        MaskedPacked& operator/=(const packed_type& rhs) noexcept {
            data = select(active_mask(), data / rhs, data);
            return *this;
        }

        MaskedPacked& operator+=(storage_type rhs) noexcept {
            data = select(active_mask(), data + rhs, data);
            return *this;
        }

        MaskedPacked& operator-=(storage_type rhs) noexcept {
            data = select(active_mask(), data - rhs, data);
            return *this;
        }

        MaskedPacked& operator*=(storage_type rhs) noexcept {
            data = select(active_mask(), data * rhs, data);
            return *this;
        }

        MaskedPacked& operator/=(storage_type rhs) noexcept {
            data = select(active_mask(), data / rhs, data);
            return *this;
        }


        // Arithmetic intentionally returns an ordinary packed value.
        friend packed_type operator+(const MaskedPacked& lhs, const MaskedPacked& rhs) noexcept {
            return lhs.data + rhs.data;
        }
        friend packed_type operator-(const MaskedPacked& lhs, const MaskedPacked& rhs) noexcept {
            return lhs.data - rhs.data;
        }
        friend packed_type operator*(const MaskedPacked& lhs, const MaskedPacked& rhs) noexcept {
            return lhs.data * rhs.data;
        }
        friend packed_type operator/(const MaskedPacked& lhs, const MaskedPacked& rhs) noexcept {
            return lhs.data / rhs.data;
        }

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
        friend packed_type operator+(const MaskedPacked& lhs, storage_type rhs) noexcept {
            return lhs.data + rhs;
        }

        friend packed_type operator+(storage_type lhs, const MaskedPacked& rhs) noexcept {
            return lhs + rhs.data;
        }

        friend packed_type operator-(const MaskedPacked& lhs, storage_type rhs) noexcept {
            return lhs.data - rhs;
        }

        friend packed_type operator-(storage_type lhs, const MaskedPacked& rhs) noexcept {
            return lhs - rhs.data;
        }

        friend packed_type operator*(const MaskedPacked& lhs, storage_type rhs) noexcept {
            return lhs.data * rhs;
        }

        friend packed_type operator*(storage_type lhs, const MaskedPacked& rhs) noexcept {
            return lhs * rhs.data;
        }

        friend packed_type operator/(const MaskedPacked& lhs, storage_type rhs) noexcept {
            return lhs.data / rhs;
        }

        friend packed_type operator/(storage_type lhs, const MaskedPacked& rhs) noexcept {
            return lhs / rhs.data;
        }


        // rotations
        template<unsigned K = 1>
        MaskedPacked& rotate_left() noexcept {
            data = data.template rotate_left<K>();
            return *this;
        }

        template<unsigned K = 1>
        MaskedPacked& rotate_right() noexcept {
            data = data.template rotate_right<K>();
            return *this;
        }

    private:
        [[nodiscard]] const mask_type& active_mask() const noexcept {
            APRIL_ASSERT(mask != nullptr, "In MaskedPacked: Mask has been bound to a valid packed mask");
            return *mask;
        }

        packed_type data = {};
        const mask_type * mask = nullptr;
    };

}
