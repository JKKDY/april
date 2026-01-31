#pragma once
#include <concepts>
#include <iterator>
#include <utility>
#include <algorithm>

namespace april::math {

    struct Range {
        struct Iterator {
            // Essential for STL algorithms
            using iterator_category = std::random_access_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = size_t;
            using pointer           = const size_t*;  // "Proxy" pointer
            using reference         = const size_t;   // Value is returned by copy usually, but const ref is ok

            size_t current;

            constexpr Iterator() : current(0) {}
            constexpr explicit Iterator(const size_t c) : current(c) {}


            // Increment/Decrement
            constexpr Iterator& operator++() noexcept {
                ++current; return *this;
            }
            constexpr Iterator operator++(int) noexcept {
                const Iterator tmp = *this; ++current; return tmp;
            }
            constexpr Iterator& operator--() noexcept {
                --current; return *this;
            }
            constexpr Iterator operator--(int) noexcept {
                const Iterator tmp = *this; --current; return tmp;
            }

            // arithmetic
            constexpr Iterator& operator+=(const difference_type n) noexcept {
                current += n; return *this;
            }
            constexpr Iterator& operator-=(const difference_type n) noexcept {
                current -= n; return *this;
            }
            constexpr difference_type operator-(const Iterator& other) const noexcept {
                return static_cast<difference_type>(current) - static_cast<difference_type>(other.current);
            }
            friend constexpr Iterator operator+(const difference_type n, const Iterator& it) noexcept {
                return Iterator{ it.current + static_cast<size_t>(n) };
            }
            friend constexpr Iterator operator+(const Iterator& it, difference_type n) noexcept {
                return Iterator{ it.current + static_cast<size_t>(n) };
            }
            friend constexpr Iterator operator-(const Iterator& it, difference_type n) noexcept {
                return Iterator{ it.current - static_cast<size_t>(n) };
            }

            // access
            constexpr size_t operator*() const noexcept { return current; }
            constexpr size_t operator[](const difference_type n) const noexcept { return current + n; }

            // ordering & equality
            constexpr auto operator<=>(const Iterator&) const noexcept = default;
            constexpr bool operator==(const Iterator&) const noexcept = default;

        };

        // Constructors
        constexpr Range() = default;

        // initialize from explicit [start, stop) range
        constexpr Range(const size_t start, const size_t stop) :
            start(start), stop(std::max(start, stop)) {}

        // initialize from a pair of integers
        template<std::integral I>
        constexpr explicit Range(const std::pair<I, I>& p)
            : Range(static_cast<size_t>(p.first), static_cast<size_t>(p.second)) {}

        // initialize from generic range (e.g. std::views::iota) using the ranges first element + size
        template<std::ranges::range R>
        requires std::convertible_to<std::ranges::range_value_t<R>, size_t>
              && std::ranges::sized_range<R>  // enforce O(1) size
              && (!std::same_as<std::remove_cvref_t<R>, Range>) // prevent eating Copy Ctor
        constexpr explicit Range(R&& r) {
            if (std::ranges::empty(r)) {
                start = stop = 0;
            } else {
                start = static_cast<size_t>(*std::ranges::begin(r));
                stop  = start + std::ranges::size(r);
            }
        }


        // Accessors & utility
        [[nodiscard]] constexpr size_t size() const noexcept {
            return stop - start;
        }
        [[nodiscard]] constexpr bool empty() const noexcept {
            return start == stop;
        }
        [[nodiscard]] constexpr bool contains(const size_t val) const noexcept {
            return val >= start && val < stop;
        }
        [[nodiscard]] constexpr bool intersects(const Range& other) const noexcept {
            return (start < other.stop) && (stop > other.start);
        }
        [[nodiscard]] constexpr Range intersection(const Range& other) const noexcept {
            const size_t s = std::max(start, other.start);
            const size_t e = std::min(stop, other.stop);
            return {s, e};
        }

        // standard iterator interface
        [[nodiscard]] constexpr Iterator begin() const noexcept { return Iterator{start}; }
        [[nodiscard]] constexpr Iterator end()   const noexcept { return Iterator{stop}; }

        // Allow array-like access to the range conceptual values
        [[nodiscard]] constexpr size_t operator[](const size_t index) const noexcept {
            return start + index;
        }

        size_t start{0};
        size_t stop{0};
    };

    // Validation
    static_assert(std::ranges::random_access_range<Range>);
    static_assert(std::ranges::sized_range<Range>);
    static_assert(std::is_trivially_copyable_v<Range>);
}