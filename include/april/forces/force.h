#pragma once

#include <variant>
#include <concepts>

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {

    struct Force {
        double cutoff;

        explicit Force(const double cutoff): cutoff(cutoff) {}

        vec3 operator()(this const auto & self,
            env::internal::Particle const& p1, env::internal::Particle const& p2, const vec3 & r) {
            static_assert(
                requires { self.eval(p1, p2, r); },
                "Force must implement eval(env::internal::Particle, env::internal::Particle, const vec3&)"
            );

            return self.eval(p1,p2,r);
        }

        auto mix_forces(this const auto& self, const auto & other) {
            static_assert(
                requires { self.mix(other); },
                "mix() not implemented"
            );

            using SelfT  = std::remove_cvref_t<decltype(self)>;
            using OtherT = std::remove_cvref_t<decltype(other)>;

            static_assert(std::same_as<SelfT, OtherT>,
                "Force::mix() requires both operands to be of the same type.");

            return self.mix(other);
        }


        [[nodiscard]] bool has_cutoff() const {
            return cutoff >= 0;
        }
    };


    // define boundary concept
    template <class F>
    concept IsForce = std::derived_from<F, Force>;

    // define Force pack
    template<IsForce... Fs>
    struct ForcePack { /* empty */ };

    template<class... Fs>
    inline constexpr ForcePack<Fs...> forces{};


    // Concept to check if a type T is a ForcePack
    template<typename T>
    inline constexpr bool is_force_pack_v = false; // Default

    template<IsForce... Fs>
    inline constexpr bool is_force_pack_v<ForcePack<Fs...>> = true; // Specialization

    template<typename T>
    concept IsForcePack = is_force_pack_v<std::remove_cvref_t<T>>;


    namespace internal {

        template<typename T>
        struct is_force_variant : std::false_type {};

        template<IsForce... Fs>
        struct is_force_variant<std::variant<Fs...>> : std::true_type {};

        template<typename T>
        concept IsForceVariant = is_force_variant<T>::value;


        template<IsForceVariant FV> struct InteractionInfo {
            bool pair_contains_types;
            std::pair<int,int> key_pair;
            FV force;

            InteractionInfo(const bool is_type_pair, const std::pair<int,int>& key, FV f)
              : pair_contains_types(is_type_pair), key_pair(key), force(std::move(f))
            {}
        };

    }
} // namespace april::env