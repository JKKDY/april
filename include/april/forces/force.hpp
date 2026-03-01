#pragma once

#include <variant>
#include <concepts>
#include <utility>

#include "april/base/types.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/scalar_access.hpp"


namespace april {
    struct NoForce;
}

namespace april::force {

    constexpr double no_cutoff = 1.0e150; // 1.0e150 squared is 1.0e300 <  max of double = 1.79e308

    enum class ForceSymmetry : uint8_t {
        Antisymmetric, // f(p1, p2) = - f(p2, p1) -> N3 applicable
        Symmetric, // f(p1, p2) = f(p2, p1) -> maybe in some esoteric active matter or graph sym
        Nonsymmetric // no relation
    };

    struct Force {
        static constexpr auto symmetry = ForceSymmetry::Antisymmetric;
        explicit Force(const double cutoff): force_cutoff(cutoff), force_cutoff2(cutoff*cutoff) {}

        // template<ParticleField IncomingMask, env::IsUserData U>
        AP_FORCE_INLINE
        auto operator()(this const auto& self,
                const auto & p1,
                const auto & p2,
                const auto & r) {

            // TODO: Force operator()
            // static_assert(
            //     requires { {self.eval(p1, p2, r)} -> std::same_as<decltype(r)>; },
            //     "Force must implement eval(auto p1, auto p2, auto r) -> decltype(r)"
            // );
            //
            // using Derived = std::remove_cvref_t<decltype(self)>;
            //
            // // check for fields requirements
            // static_assert(
            //     requires { Derived::fields; },
            //     "Force subclass must define 'static constexpr env::Field fields'"
            // );
            //
            // constexpr ParticleField Required = Derived::fields;
            //
            // // check for
            // static_assert(
            //     (IncomingMask & Required) == Required,
            //     "ParticleView is missing required fields for this Force."
            // );

            return self.eval(p1, p2, r);
        }

        auto mix_forces(this const auto& self, const auto & other) {
            using SelfT  = std::remove_cvref_t<decltype(self)>;
            using OtherT = std::remove_cvref_t<decltype(other)>;

            static_assert(std::same_as<SelfT, OtherT>,
                "Force::mix() requires both operands to be of the same type.");

            if constexpr (requires{self.mix(other);}) {
                return self.mix(other);
            } else {
                if (!self.equals(other)) {
                    throw std::invalid_argument("Mixing disabled by default for this force");
                } else {
                    return self;
                }
            }
        }

        [[nodiscard]] bool has_cutoff() const noexcept{
            return cutoff() < no_cutoff;
        }

        [[nodiscard]] double cutoff() const noexcept{
            return force_cutoff;
        }
        [[nodiscard]] double cutoff2() const noexcept{
            return force_cutoff2;
        }

        auto with_cutoff(this auto && self, const double c) {
            self.force_cutoff = c;
            self.force_cutoff2 = c*c;
            return self;
        }

        bool equals(this const auto & self, const auto & other) {
            using SelfT  = std::remove_cvref_t<decltype(self)>;
            using OtherT = std::remove_cvref_t<decltype(other)>;

            if constexpr (!std::same_as<SelfT, OtherT>) {
                return false;
            }

            if (self.force_cutoff != other.force_cutoff) {
                return false;
            }

            return self == other;
        }

        bool operator==(const Force&) const = default;

    protected:
        double force_cutoff;
        double force_cutoff2;
    };


    // define boundary concept
    template <class F>
    concept IsForce = std::derived_from<F, Force>;


    namespace internal {
        // define Force pack
        template<IsForce... Fs> struct ForcePack {};

        // Concept to check if a type T is a ForcePack
        template<typename T>
        inline constexpr bool is_force_pack_v = false; // Default

        template<IsForce... Fs>
        inline constexpr bool is_force_pack_v<ForcePack<Fs...>> = true; // Specialization

        template<typename T>
        concept IsForcePack = is_force_pack_v<std::remove_cvref_t<T>>;


        // check if std::variant of forces
        template<typename T>
        struct is_force_variant : std::false_type {};

        template<IsForce... Fs>
        struct is_force_variant<std::variant<Fs...>> : std::true_type {};

        template<typename T>
        concept IsForceVariant = is_force_variant<T>::value;


        template<IsForceVariant FV> struct TypeInteraction {
            const ParticleType type1;
            const ParticleType type2;
            const FV force;

            TypeInteraction(const ParticleType type1, const ParticleType type2, FV f)
              : type1(std::min(type1, type2)), type2(std::max(type1, type2)), force(std::move(f))
            {}
        };

        template<IsForceVariant FV> struct IdInteraction {
            const ParticleID id1;
            const ParticleID id2;
            const FV force;

            IdInteraction(const ParticleID id1, const ParticleID id2, FV f)
              : id1(std::min(id1, id2)), id2(std::max(id1, id2)), force(std::move(f))
            {}
        };


        // internal placeholder only
        struct ForceSentinel : Force {
            static constexpr auto fields = ParticleField::none;

            ForceSentinel() : Force(-1.0) {}

            template<ParticleField M, particle::IsParticleAttributes U>
            vec3 eval(const particle::internal::ScalarParticleView<M, U> &, const particle::internal::ScalarParticleView<M, U> &, const vec3&) const noexcept {
                AP_ASSERT(false, "NullForce should never be executed");
                std::unreachable();
            }
            [[nodiscard]] ForceSentinel mix(ForceSentinel const&) const { return {}; }

            bool operator==(const ForceSentinel&) const = default;
        };


        template<class... Fs>
        struct VariantType {
            // Disallow the internal sentinel type in user packs
            static_assert((!std::is_same_v<ForceSentinel, Fs> && ...),
                          "ForceSentinel must NOT appear in ForcePack (internal sentinel only).");

            // Detect whether NoForce is already supplied
            static constexpr bool has_no_force = (std::is_same_v<NoForce, Fs> || ...);

            // Compute the variant type
            using type = std::conditional_t<
                has_no_force,
                std::variant<ForceSentinel, Fs...>,           // user already included NoForce
                std::variant<ForceSentinel, Fs..., NoForce>   // append it
            >;
        };

        // Convenience alias
        template<class... Fs>
        using VariantType_t = VariantType<Fs...>::type;
    } // namespace internal
} // namespace april::force


namespace april {
    // define Force pack
    template<class... Fs> inline constexpr force::internal::ForcePack<Fs...> forces{};
}









