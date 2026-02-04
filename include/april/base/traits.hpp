#pragma once

#include <concepts>

namespace april {
    template <typename T, typename... Ts>
    concept same_as_any = (... or std::same_as<T, Ts>);

    template <typename T>
    concept IsScalar =  std::floating_point<T> || std::integral<T>;

}