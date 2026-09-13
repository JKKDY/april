#pragma once

#include <concepts>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <algorithm>


#include "april/base/types.hpp"
#include "april/particle/access/packed_access.hpp"
#include "april/particle/properties.hpp"
#include "april/exec/policy.hpp"

namespace april {
    struct NoInteraction;
}

namespace april::interaction {

    constexpr double no_cutoff = 1.0e150; // 1.0e150 squared is 1.0e300 <  max of double = 1.79e308

    enum class InteractionSymmetry : uint8_t {
        Antisymmetric, // f(p1, p2) = - f(p2, p1) -> N3 applicable
        Symmetric, // f(p1, p2) = f(p2, p1) -> maybe in some esoteric active matter or graph sim
        Nonsymmetric // no relation
    };

    /**
     * @brief Base class for pairwise particle interactions.
     *
     * Custom interactions derive from Interaction, declare the particle fields they need through
     * `static constexpr ParticleField fields`, and implement `eval(p1, p2, r)`.
     * The same implementation may be used for scalar and SIMD particle views.
     */
    struct Interaction {
        static constexpr auto symmetry = InteractionSymmetry::Antisymmetric;
        static constexpr auto vector_mode = exec::ExecutionMode::Scalar | exec::ExecutionMode::Packed; // scalar only must be a deliberate opt-out

        explicit Interaction(const double cutoff): interaction_cutoff(cutoff), interaction_cutoff2(cutoff*cutoff) {}

        APRIL_FORCE_INLINE
        auto operator()(this const auto& self, const auto & p1, const auto & p2, const auto & r) {
            using Derived = std::remove_cvref_t<decltype(self)>;
            using ReturnType = std::remove_cvref_t<decltype(r)>; // Resolves to vec3 or pvec3

            // check if fields is defined
            static_assert(
                 requires {
                     []{ constexpr auto _ = Derived::fields; };
                     { Derived::fields } -> std::convertible_to<ParticleField>;
                 },
                 "[APRIL] Interaction: subclass must define 'static constexpr ParticleField fields'"
             );

            // check if the incoming particles have all requested fields
            using P1Type = std::remove_cvref_t<decltype(p1)>;
            constexpr ParticleField Required = Derived::fields;
            constexpr ParticleField IncomingMask = P1Type::ReadAccess;
            static_assert(
                (IncomingMask & Required) == Required,
                "[APRIL] Interaction: ParticleView is missing required fields for this Interaction."
            );


            constexpr bool is_vector = particle::IsPackedParticleAccessor<P1Type>;
            constexpr auto packed_mode = exec::ExecutionMode::Packed;
            constexpr auto scalar_mode = exec::ExecutionMode::Scalar;
            if constexpr (is_vector) {
                // Try Vector Override First
                if constexpr (requires { self.eval_vector(p1, p2, r); }) {
                    static_assert(requires { { self.eval_vector(p1, p2, r) } -> std::same_as<ReturnType>; },
                        "[APRIL] Interaction: eval_vector must return the same vector type as 'r' (pvec3)");
                    return self.eval_vector(p1, p2, r);
                }
                // Try Templated Generic Next
                else if constexpr (requires { self.template eval<packed_mode>(p1, p2, r); }) {
                    static_assert(requires { { self.template eval<packed_mode>(p1, p2, r) } -> std::same_as<ReturnType>; },
                        "[APRIL] Interaction: eval<true> must return the same vector type as 'r' (pvec3)");
                    return self.template eval<packed_mode>(p1, p2, r);
                }
                // Fallback to Non-Templated Generic
                else {
                    static_assert(requires { { self.eval(p1, p2, r) } -> std::same_as<ReturnType>; },
                        "[APRIL] Interaction: must implement eval(p1, p2, r), eval<packed>(p1, p2, r), or eval_vector(p1, p2, r)");
                    return self.eval(p1, p2, r);
                }
            } else {
                // Try Scalar Override First
                if constexpr (requires { self.eval_scalar(p1, p2, r); }) {
                    static_assert(requires { { self.eval_scalar(p1, p2, r) } -> std::same_as<ReturnType>; },
                        "[APRIL] Interaction: eval_scalar must return the same vector type as 'r' (vec3)");
                    return self.eval_scalar(p1, p2, r);
                }
                // Try Templated Generic Next
                else if constexpr (requires { self.template eval<scalar_mode>(p1, p2, r); }) {
                    static_assert(requires { { self.template eval<scalar_mode>(p1, p2, r) } -> std::same_as<ReturnType>; },
                        "[APRIL] Interaction: eval<false> must return the same vector type as 'r' (vec3)");
                    return self.template eval<scalar_mode>(p1, p2, r);
                }
                // Fallback to Non-Templated Generic
                else {
                    static_assert(requires { { self.eval(p1, p2, r) } -> std::same_as<ReturnType>; },
                        "[APRIL] Interaction: must implement eval(p1, p2, r), eval<packed>(p1, p2, r), or eval_scalar(p1, p2, r)");
                    return self.eval(p1, p2, r);
                }
            }
        }

        auto mix_interactions(this const auto& self, const auto & other) {
            using SelfT  = std::remove_cvref_t<decltype(self)>;
            using OtherT = std::remove_cvref_t<decltype(other)>;

            static_assert(std::same_as<SelfT, OtherT>,
                "[APRIL] Error: Interaction::mix_interactions() requires both operands to be of the same type.");

            if constexpr (requires{self.mix(other);}) {
                return self.mix(other);
            } else {
                if (!self.equals(other)) {
                    throw std::invalid_argument("[APRIL] Error: Mixing disabled by default for this interaction");
                } else {
                    return self;
                }
            }
        }

        [[nodiscard]] bool has_cutoff() const noexcept{
            return cutoff() < no_cutoff;
        }

        [[nodiscard]] double cutoff() const noexcept{
            return interaction_cutoff;
        }
        [[nodiscard]] double cutoff2() const noexcept{
            return interaction_cutoff2;
        }

        auto&& with_cutoff(this auto && self, const double c) {
            self.interaction_cutoff = c;
            self.interaction_cutoff2 = c*c;
            return self;
        }

        bool equals(this const auto & self, const auto & other) {
            using SelfT  = std::remove_cvref_t<decltype(self)>;
            using OtherT = std::remove_cvref_t<decltype(other)>;

            if constexpr (!std::same_as<SelfT, OtherT>) {
                return false;
            }

            if (self.interaction_cutoff != other.interaction_cutoff) {
                return false;
            }

            return self == other;
        }

        bool operator==(const Interaction&) const = default;

    protected:
        double interaction_cutoff;
        double interaction_cutoff2;
    };


    // Define the pairwise interaction concept.
    template <class I>
    concept IsInteraction = std::derived_from<I, Interaction>;


    namespace internal {
        // Define the interaction pack used by the public interactions<T...> marker.
        template<IsInteraction... InteractionTs> struct InteractionPack {};

        // Concept to check if a type T is a InteractionPack
        template<typename T>
        inline constexpr bool is_interaction_pack_v = false; // Default

        template<IsInteraction... InteractionTs>
        inline constexpr bool is_interaction_pack_v<InteractionPack<InteractionTs...>> = true; // Specialization

        template<typename T>
        concept IsInteractionPack = is_interaction_pack_v<std::remove_cvref_t<T>>;


        // Check if std::variant contains only interactions.
        template<typename T>
        struct is_interaction_variant : std::false_type {};

        template<IsInteraction... InteractionTs>
        struct is_interaction_variant<std::variant<InteractionTs...>> : std::true_type {};

        template<typename T>
        concept IsInteractionVariant = is_interaction_variant<T>::value;


        template<IsInteractionVariant IV> struct TypeInteraction {
            const ParticleType type1;
            const ParticleType type2;
            const IV interaction;

            TypeInteraction(const ParticleType type1, const ParticleType type2, IV value)
              : type1(std::min(type1, type2)), type2(std::max(type1, type2)), interaction(std::move(value))
            {}
        };

        template<IsInteractionVariant IV> struct IdInteraction {
            const ParticleID id1;
            const ParticleID id2;
            const IV interaction;

            IdInteraction(const ParticleID id1, const ParticleID id2, IV value)
              : id1(std::min(id1, id2)), id2(std::max(id1, id2)), interaction(std::move(value))
            {}
        };


        // internal placeholder only
        struct InteractionSentinel : Interaction {
            static constexpr auto fields = ParticleField::none;

            InteractionSentinel() : Interaction(-1.0) {}

            vec3 eval(auto, auto, const vec3&) const noexcept {
                APRIL_ASSERT(false, "InteractionSentinel should never be executed");
                std::unreachable();
            }
            [[nodiscard]] InteractionSentinel mix(InteractionSentinel const&) const { return {}; }

            bool operator==(const InteractionSentinel&) const = default;
        };


        template<class... InteractionTs>
        struct InteractionVariant {
            // Disallow the internal sentinel type in user packs
            static_assert((!std::is_same_v<InteractionSentinel, InteractionTs> && ...),
                          "[APRIL] Error: InteractionSentinel must NOT appear in InteractionPack (internal sentinel only).");

            // Detect whether NoInteraction is already supplied
            static constexpr bool has_no_interaction = (std::is_same_v<NoInteraction, InteractionTs> || ...);

            // Compute the variant type
            using type = std::conditional_t<
                has_no_interaction,
                std::variant<InteractionSentinel, InteractionTs...>,
                std::variant<InteractionSentinel, InteractionTs..., NoInteraction>
            >;
        };

        // Convenience alias
        template<class... InteractionTs>
        using InteractionVariant_t = InteractionVariant<InteractionTs...>::type;
    } // namespace internal
} // namespace april::interaction


namespace april {
    // Public interaction declaration marker.
    template<class... InteractionTs>
    inline constexpr interaction::internal::InteractionPack<InteractionTs...> interactions{};
}











