#pragma once

#include <variant>
#include <concepts>
#include <utility>
#include <limits>

#include "april/common.h"
#include "april/particle/fields.h"


namespace april::force {

    constexpr double no_cutoff = std::numeric_limits<double>::max();

    struct Force {
        // TODO: should this/can this be const?
        double cutoff;
        double cutoff2;

        explicit Force(const double cutoff): cutoff(cutoff), cutoff2(cutoff*cutoff) {}

        template<env::IsConstFetcher F>
        vec3 operator()(this const auto & self, const F & p1, const F & p2, const vec3 & r) {
            static_assert(
                requires { {self.eval(p1, p2, r)} -> std::same_as<vec3>; },
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
            return cutoff != no_cutoff;
        }
    };


    // define boundary concept
    template <class F>
    concept IsForce = std::derived_from<F, Force>;

    // define Force pack
    template<IsForce... Fs> struct ForcePack {};
    template<class... Fs> inline constexpr ForcePack<Fs...> forces{};


    // Concept to check if a type T is a ForcePack
    template<typename T>
    inline constexpr bool is_force_pack_v = false; // Default

    template<IsForce... Fs>
    inline constexpr bool is_force_pack_v<ForcePack<Fs...>> = true; // Specialization

    template<typename T>
    concept IsForcePack = is_force_pack_v<std::remove_cvref_t<T>>;


    struct NoForce;

    namespace internal {

        template<typename T>
        struct is_force_variant : std::false_type {};

        template<IsForce... Fs>
        struct is_force_variant<std::variant<Fs...>> : std::true_type {};

        template<typename T>
        concept IsForceVariant = is_force_variant<T>::value;


        template<IsForceVariant FV> struct TypeInteraction {
            const env::ParticleType type1;
            const env::ParticleType type2;
            const FV force;

            TypeInteraction(const env::ParticleType type1, const env::ParticleType type2, FV f)
              : type1(std::min(type1, type2)), type2(std::max(type1, type2)), force(std::move(f))
            {}
        };

        template<IsForceVariant FV> struct IdInteraction {
            const env::ParticleID id1;
            const env::ParticleID id2;
            const FV force;

            IdInteraction(const env::ParticleID id1, const env::ParticleID id2, FV f)
              : id1(std::min(id1, id2)), id2(std::max(id1, id2)), force(std::move(f))
            {}
        };


        // internal placeholder only
        struct ForceSentinel : Force {
            ForceSentinel() : Force(-1.0) {}

		    template<env::IsConstFetcher F>
            vec3 eval(const F&, const F&, const vec3&) const noexcept {
                AP_ASSERT(false, "NullForce should never be executed");
                std::unreachable();
            }
            [[nodiscard]] ForceSentinel mix(ForceSentinel const&) const { return {}; }
        };


        template<class... Fs>
        struct VariantType {
            // 1. Disallow the internal sentinel type in user packs
            static_assert((!std::is_same_v<ForceSentinel, Fs> && ...),
                          "ForceSentinel must NOT appear in ForcePack (internal sentinel only).");

            // 2. Detect whether NoForce is already supplied
            static constexpr bool has_no_force = (std::is_same_v<NoForce, Fs> || ...);

            // 3. Compute the variant type
            using type = std::conditional_t<
                has_no_force,
                std::variant<ForceSentinel, Fs...>,           // user already included NoForce
                std::variant<ForceSentinel, Fs..., NoForce>   // append it
            >;
        };

        // Convenience alias
        template<class... Fs>
        using VariantType_t = typename VariantType<Fs...>::type;
    }
} // namespace april::env