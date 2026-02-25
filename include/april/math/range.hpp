#pragma once
#include <concepts>
#include <iterator>
#include <utility>
#include <algorithm>
#include <type_traits>

namespace april::math {

    struct Range {
        struct Iterator {
            // C++20 compliant tags for iterators returning by value
            using iterator_concept  = std::random_access_iterator_tag;
            using iterator_category = std::input_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = size_t;
            using pointer           = void;
            using reference         = size_t;

            size_t current;
            size_t step;

            constexpr Iterator() : current(0), step(1) {}
            constexpr explicit Iterator(const size_t c, const size_t s) : current(c), step(s) {}

            // Increment/Decrement
            constexpr Iterator& operator++() noexcept {
                current += step; return *this;
            }
            constexpr Iterator operator++(int) noexcept {
                const Iterator tmp = *this; current += step; return tmp;
            }
            constexpr Iterator& operator--() noexcept {
                current -= step; return *this;
            }
            constexpr Iterator operator--(int) noexcept {
                const Iterator tmp = *this; current -= step; return tmp;
            }

            // Safe Arithmetic
            constexpr Iterator& operator+=(const difference_type n) noexcept {
                if (n >= 0) current += static_cast<size_t>(n) * step;
                else        current -= static_cast<size_t>(-n) * step;
                return *this;
            }
            constexpr Iterator& operator-=(const difference_type n) noexcept {
                return *this += -n;
            }
            constexpr difference_type operator-(const Iterator& other) const noexcept {
                if (current >= other.current) {
                    return static_cast<difference_type>((current - other.current) / step);
                } else {
                    return -static_cast<difference_type>((other.current - current) / step);
                }
            }
            friend constexpr Iterator operator+(const difference_type n, Iterator it) noexcept {
                it += n; return it;
            }
            friend constexpr Iterator operator+(Iterator it, const difference_type n) noexcept {
                it += n; return it;
            }
            friend constexpr Iterator operator-(Iterator it, const difference_type n) noexcept {
                it -= n; return it;
            }

            // Access
            constexpr size_t operator*() const noexcept { return current; }
            constexpr size_t operator[](const difference_type n) const noexcept {
                return *(*this + n);
            }

            // Ordering & Equality
            constexpr auto operator<=>(const Iterator&) const noexcept = default;
            constexpr bool operator==(const Iterator&) const noexcept = default;
        };

        // Constructors
        constexpr Range() = default;

        // Initialize from [start, stop) range with optional step
        constexpr Range(const size_t start, const size_t stop, const size_t step = 1) :
            start(start), stop(std::max(start, stop)), step(step > 0 ? step : 1) {}

        // Initialize from a pair of integers
        template<std::integral I>
        constexpr explicit Range(const std::pair<I, I>& p)
            : Range(static_cast<size_t>(p.first), static_cast<size_t>(p.second)) {}

        // Accessors & Utility
        [[nodiscard]] constexpr size_t size() const noexcept {
            if (stop <= start) return 0;
            return (stop - start - 1) / step + 1;
        }
        [[nodiscard]] constexpr bool empty() const noexcept {
            return size() == 0;
        }
        [[nodiscard]] constexpr bool contains(const size_t val) const noexcept {
            return val >= start && val < stop && ((val - start) % step == 0);
        }

        // Standard iterator interface
        [[nodiscard]] constexpr Iterator begin() const noexcept { return Iterator{start, step}; }
        [[nodiscard]] constexpr Iterator end()   const noexcept {
            return Iterator{start + size() * step, step};
        }

        // Allow array-like access to the range conceptual values
        [[nodiscard]] constexpr size_t operator[](const size_t index) const noexcept {
            return start + index * step;
        }

        size_t start{0};
        size_t stop{0};
        size_t step{1};
    };

    // Validation
    static_assert(std::ranges::random_access_range<Range>);
    static_assert(std::ranges::sized_range<Range>);
    static_assert(std::is_trivially_copyable_v<Range>);
}