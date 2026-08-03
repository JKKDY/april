#pragma once

#include <concepts>

namespace april {
    template <typename T, typename... Ts>
    concept same_as_any = (... or std::same_as<T, Ts>);

    template<typename T>
    concept IsScalar =
    !std::same_as<std::remove_cvref_t<T>, bool> &&
    (std::floating_point<std::remove_cvref_t<T>> ||
     std::integral<std::remove_cvref_t<T>>);
}














