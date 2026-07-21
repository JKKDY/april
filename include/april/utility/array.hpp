#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <ranges>

namespace april::utility {
    template <auto A, auto B>
    consteval auto array_intersection()
    {
        using T = std::ranges::range_value_t<decltype(A)>;
        static_assert(std::same_as<T, std::ranges::range_value_t<decltype(B)>>);
        static_assert(std::equality_comparable<T>);

        constexpr std::size_t capacity = std::min(A.size(), B.size());

        struct scratch_result {
            std::array<T, capacity> values{};
            std::size_t count{};
        };

        constexpr auto scratch = [] {
            scratch_result result{};

            for (T value : A) {
                const auto used = result.values | std::views::take(result.count);

                if (std::ranges::contains(B, value) &&
                    !std::ranges::contains(used, value)) {
                    result.values[result.count++] = value;
                    }
            }

            return result;
        }();

        std::array<T, scratch.count> exact{};

        std::ranges::copy(
            scratch.values | std::views::take(scratch.count),
            exact.begin()
        );

        return exact;
    }
}
