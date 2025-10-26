#pragma once

#include <variant>
#include <concepts>

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {

    // TODO define a force superclass instead
    template<typename F> concept IsForce =
        std::copy_constructible<F> &&
        std::assignable_from<F&, F const&> &&
        std::movable<F> &&
    requires(F const& f,
             env::internal::Particle const& p1,
             env::internal::Particle const& p2,
             vec3 const& r,
             F const& o)
    {
        { f(p1, p2, r) } noexcept -> std::same_as<vec3>;
        { f.mix(o) } -> std::same_as<F>;
        { f.cutoff_radius } -> std::convertible_to<double>;
    };


    template<IsForce... Fs>
    struct ForcePack { /* empty */ };

    template<class... Fs>
    inline constexpr ForcePack<Fs...> forces{};


    namespace internal {

        template<typename T>
        struct is_force_variant : std::false_type {};

        template<IsForce... Fs>
        struct is_force_variant<std::variant<Fs...>> : std::true_type {};

        template<typename T>
        concept ForceVariant = is_force_variant<T>::value;


        template<ForceVariant FV> struct InteractionInfo {
            bool pair_contains_types;
            std::pair<int,int> key_pair;
            FV force;

            InteractionInfo(const bool is_type_pair, const std::pair<int,int>& key, FV f)
              : pair_contains_types(is_type_pair), key_pair(key), force(std::move(f))
            {}
        };

    }
} // namespace april::env