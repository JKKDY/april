#pragma once

#include <variant>
#include <concepts>
#include <utility>

#include "april/base/types.hpp"
#include "april/particle/fields.hpp"
#include "april/particle/defs.hpp"
#include "april/particle/access.hpp"


namespace april::force {

    constexpr double no_cutoff = 1.0e150; // 1.0e150 squared is 1.0e300 <  max of double = 1.79e308

    struct Force {
        explicit Force(const double cutoff): force_cutoff(cutoff), force_cutoff2(cutoff*cutoff) {}

        template<env::FieldMask IncomingMask, env::IsUserData U>
        vec3 operator()(this const auto& self,
                const env::ParticleView<IncomingMask, U> & p1,
                const env::ParticleView<IncomingMask, U> & p2,
                const vec3 & r) {

            static_assert(
                requires { {self.eval(p1, p2, r)} -> std::same_as<vec3>; },
                "Force must implement eval(env::internal::Particle, env::internal::Particle, const vec3&)"
            );


            using Derived = std::remove_cvref_t<decltype(self)>;

            // check for fields requirements
            static_assert(
                requires { Derived::fields; },
                "Force subclass must define 'static constexpr env::FieldMask fields'"
            );

            constexpr env::FieldMask Required = Derived::fields;

            // check for
            static_assert(
                (IncomingMask & Required) == Required,
                "ParticleView is missing required fields for this Force."
            );

            return self.eval(p1, p2, r);
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
            static constexpr env::FieldMask fields = +env::Field::none;

            ForceSentinel() : Force(-1.0) {}

            template<env::FieldMask M, env::IsUserData U>
            vec3 eval(const env::ParticleView<M, U> &, const env::ParticleView<M, U> &, const vec3&) const noexcept {
                AP_ASSERT(false, "NullForce should never be executed");
                std::unreachable();
            }
            [[nodiscard]] ForceSentinel mix(ForceSentinel const&) const { return {}; }

            bool operator==(const ForceSentinel&) const = default;
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