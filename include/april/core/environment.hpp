/**
 * @file environment.hpp
 * @brief User-facing staging area for simulation configuration.
 *
 * An Environment collects particles, interactions, fields, boundaries,
 * controllers, and domain information before they are materialized into
 * a simulation-ready System by build_system(...).
 */

#pragma once

#include <any>
#include <array>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "april/base/types.hpp"
#include "april/core/domain.hpp"
#include "april/core/internal/environment_traits.hpp"

#include "april/interactions/force.hpp"
#include "april/boundaries/boundary.hpp"
#include "april/controllers/controller.hpp"
#include "april/fields/field.hpp"

#include "april/particle/generators.hpp"
#include "april/particle/properties.hpp"
#include "april/particle/attributes.hpp"


namespace april {

    // scoping markers for interactions
    /// @brief Target a single particle type for self-interaction.
    struct to_type {
        ParticleType type;
        constexpr explicit to_type(ParticleType t) noexcept : type(t) {}
    };

    /// @brief Target interactions between two particle types.
    struct between_types {
        ParticleType t1, t2;
        constexpr between_types(ParticleType a, ParticleType b) noexcept : t1(a), t2(b) {}
    };

    /// @brief Target interactions between two specific particles.
    struct between_ids {
        ParticleID id1, id2;
        constexpr between_ids(ParticleID a, ParticleID b) noexcept : id1(a), id2(b) {}
    };

    /**
      * @brief Declarative staging object for simulation setup.
      *
      * Environment stores the particles, interactions, fields, boundaries,
      * controllers, and domain description that are later materialized into a
      * simulation-ready System by build_system(...).
      *
      * @tparam FPack Forces pack: the set of force types supported by this environment.
      * @tparam BPack Boundaries pack: the set of boundary types supported by this environment.
      * @tparam CPack Controllers pack: the set of controller types supported by this environment.
      * @tparam FFPack Global Fields pack: the set of force field types supported by this environment.
      * @tparam ParticleData User-defined particle attributes.
      */
    template<
        interactions::internal::IsForcePack FPack,
        boundary::internal::IsBoundaryPack BPack,
        controller::internal::IsControllerPack CPack,
        field::internal::IsFieldPack FFPack,
        particle::IsParticleAttributes ParticleData>
    class Environment {
    public:
        using traits = core::internal::EnvironmentTraits<FPack, BPack, CPack, FFPack, ParticleData>;

        explicit Environment(FPack, BPack, CPack, FFPack, ParticleData) {}
        Environment() : Environment(forces<>, boundaries<>, controllers<>, fields<>, NoParticleAttributes{}) {}

        /**
         * @brief Deduce Environment template parameters from component type packs.
         *
         * Enables:
         * auto env = Environment(forces<LennardJones>, boundaries<ReflectiveBoundary>);
         */
        template<class... Args>
        requires (core::internal::is_any_pack_v<std::remove_cvref_t<Args>> && ...) &&
            (!std::same_as<std::remove_cvref_t<Args>, Environment> && ...)
        explicit Environment(Args&&...)
            : Environment(
                core::internal::get_pack_t<interactions::internal::ForcePack, Args...>{},
                core::internal::get_pack_t<boundary::internal::BoundaryPack, Args...>{},
                core::internal::get_pack_t<controller::internal::ControllerPack,Args...>{},
                core::internal::get_pack_t<field::internal::FieldPack, Args...>{},
                core::internal::get_particle_attributes_t<Args...>{}
            ) {}

    private:
        traits::environment_data_t data;
        friend auto core::internal::get_env_data<FPack, BPack, CPack, FFPack> (const Environment& env);

    public:

        //--------------
        // ADD PARTICLES
        //--------------
        /// @brief Add a single particle to the environment.
        void add_particle(const Particle& particle) {
            if (particle.id.has_value() && data.user_particle_ids.contains(particle.id.value())) {
                throw std::invalid_argument("specified id is not unique");
            }
            data.particles.push_back(particle);

            // register particle type and id (if available)
            data.user_particle_types.insert(particle.type);
            if (particle.id.has_value()) {
                data.user_particle_ids.insert(particle.id.value());
            }
        }

        /// @brief Construct and add a particle from raw scalars.
        void add_particle(
            const vec3& position,
            const vec3& velocity,
            const double mass,
            const ParticleType type = 0,
            const std::optional<ParticleID> id = {},
            const std::any& user_data = {}
        ) {
            Particle p;
            p.id = id;
            p.type = type;
            p.position = position;
            p.velocity = velocity;
            p.mass = mass;
            p.state = ParticleState::ALIVE;
            p.user_data = std::move(user_data);
            add_particle(p);
        }

        /// @brief Add a collection of particles.
        void add_particles(const std::vector<Particle>& particles) {
            this->data.particles.reserve(this->data.particles.size() + particles.size());
            for (auto & p : particles)  add_particle(p);
        }

        /// @brief Add particles generated by a primitive (e.g., Cuboid, Sphere).
        template<IsParticleGenerator T>
        void add_particles(const T& generator) {
            add_particles(generator.to_particles());
        }


        //-----------------
        // ADD INTERACTIONS
        //-----------------
        /// @brief Define self-interaction for a specific particle type.
        template<interactions::IsForce F> requires traits::template is_valid_force_v<F>
        void add_interaction(F force, to_type scope) {
            data.type_interactions.emplace_back(scope.type, scope.type, typename traits::force_variant_t{std::move(force)});
        }

        /// @brief Define interaction between two distinct particle types.
        template<interactions::IsForce F> requires traits::template is_valid_force_v<F>
        void add_interaction(F force, between_types scope) {
            data.type_interactions.emplace_back(scope.t1, scope.t2, typename traits::force_variant_t{std::move(force)});
        }

        /// @brief Define interaction between two specific particle IDs (Bonds).
        template<interactions::IsForce F> requires traits::template is_valid_force_v<F>
        void add_interaction(F force, between_ids scope) {
            data.id_interactions.emplace_back(scope.id1, scope.id2, typename traits::force_variant_t{std::move(force)});
        }


        //---------------
        // ADD BOUNDARIES
        //---------------
        /// @brief Apply a boundary condition to a specific face of the domain.
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundary(B boundary, const DomainFace face) {
            data.boundaries[boundary::face_to_int(face)].template emplace<B>(std::move(boundary));
        }

        /// @brief Apply a boundary condition to multiple faces simultaneously.
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundaries(B boundary, const std::vector<DomainFace> & faces) {
            for (const DomainFace face : faces) {
                data.boundaries[boundary::face_to_int(face)].template emplace<B>(boundary);
            }
        }

        /// @brief Apply per-face boundaries from an array. Order is -X, +X, -Y, +Y, -Z, +Z
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundaries(const std::array<B, 6> & boundaries) {
            for (const DomainFace face : all_faces) {
                data.boundaries[boundary::face_to_int(face)].template emplace<B>(boundaries[boundary::face_to_int(face)]);
            }
        }


        //----------------
        // ADD CONTROLLERS
        //----------------
        /// @brief Add a controller (e.g., Thermostat) to the system.
        template<controller::IsController C> requires traits::template is_valid_controller_v<std::remove_cvref_t<C>>
        void add_controller(C && controller) {
            data.controllers.add(controller);
        }

        /// @brief Add multiple controllers
        template<controller::IsController... C>
        void add_controllers(C&&... controllers) {
            (data.controllers.add(std::forward<C>(controllers)), ...);
        }



        //-----------
        // ADD FIELDS
        //-----------
        /// @brief Add a global force field (e.g., Uniform Gravity).
        template<field::IsField F>  requires traits::template is_valid_field_v<std::remove_cvref_t<F>>
        void add_field(F && field) {
            data.fields.add(field);
        }

        /// @brief Add multiple fields
        template<field::IsField... F>
        void add_fields(F&&... fields) {
            (data.fields.add(std::forward<F>(fields)), ...);
        }


        //-----------
        // SET DOMAIN
        //-----------
        /// @brief Set the coordinate of the minimum corner.
        void set_origin(const vec3d& origin) {
            this->data.domain.origin = origin;
        }
        /// @brief Set the coordinate of the minimum corner from scalars.
        void set_origin(const double x, const double y, const double z) {
            set_origin({x,y,z});
        }

        /// @brief Set the dimensions of the simulation box.
        void set_extent(const vec3d& extent) {
            this->data.domain.extent = extent;
        }
        /// @brief Set the dimensions of the simulation box from scalars.
        void set_extent(const double x, const double y, const double z) {
            set_extent({x,y,z});
        }

        /// @brief Set the full domain geometry.
        void set_domain(const Domain& domain) {
            data.domain = domain;
        }

        /**
         * @brief Define a fixed minimum distance between the particles and the domain boundaries.
         * @details During the build phase, the domain will be at least the particle AABB
         * expanded by this absolute value in all directions.
         */
        void domain_padding(const vec3d& margin_abs) {
            data.margin_abs = margin_abs;
        }
        void domain_padding(const double margin_abs) {
            data.margin_abs = vec3d{margin_abs};
        }

        /**
         * @brief Define a relative margin based on the current particle distribution size.
         * @details For example, a factor of 0.5 adds 50% of the particle distribution's
         * width as padding to each side.
         * @note If both margin_abs and margin_fac are set, the larger of the two is used per axis.
         */
        void domain_padding_factor(const vec3d& margin_fac) {
            data.margin_fac = margin_fac;
        }

        void domain_padding_factor(const double margin_fac) {
            domain_padding_factor(vec3d{margin_fac});
        }


        //-------------------
        // DSL-STYLE CHAINING
        //-------------------
        auto&& with_particle(this auto&& self, const Particle& p) {
            self.add_particle(p);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_particle(
            this auto&& self,
            const vec3& position,
            const vec3& velocity,
            double mass,
            ParticleType type = 0,
            std::optional<ParticleID> id = {},
            std::any user_data = {}
        ) {
            self.add_particle(position, velocity, mass, type, id, std::move(user_data));
            return std::forward<decltype(self)>(self);
        }

        auto&& with_particles(this auto&& self, const std::vector<Particle>& ps) {
            self.add_particles(ps);
            return std::forward<decltype(self)>(self);
        }

        template<IsParticleGenerator G>
        auto&& with_particles(this auto&& self, const G& particles) {
            self.add_particles(particles);
            return std::forward<decltype(self)>(self);
        }

        template<interactions::IsForce F>
            requires traits::template is_valid_force_v<std::remove_cvref_t<F>>
        auto&& with_interaction(this auto&& self, F&& force, to_type scope) {
            self.add_interaction(std::forward<F>(force), scope);
            return std::forward<decltype(self)>(self);
        }

        template<interactions::IsForce F>
            requires traits::template is_valid_force_v<std::remove_cvref_t<F>>
        auto&& with_interaction(this auto&& self, F&& force, between_types scope) {
            self.add_interaction(std::forward<F>(force), scope);
            return std::forward<decltype(self)>(self);
        }

        template<interactions::IsForce F>
            requires traits::template is_valid_force_v<std::remove_cvref_t<F>>
        auto&& with_interaction(this auto&& self, F&& force, between_ids scope) {
            self.add_interaction(std::forward<F>(force), scope);
            return std::forward<decltype(self)>(self);
        }

        template<boundary::IsBoundary B>
            requires traits::template is_valid_boundary_v<std::remove_cvref_t<B>>
        auto&& with_boundary(this auto&& self, B&& boundary, DomainFace face) {
            self.set_boundary(std::forward<B>(boundary), face);
            return std::forward<decltype(self)>(self);
        }

        template<boundary::IsBoundary B>
            requires traits::template is_valid_boundary_v<std::remove_cvref_t<B>>
        auto&& with_boundaries(this auto&& self, B&& boundary, const std::vector<DomainFace>& faces) {
            self.set_boundaries(std::forward<B>(boundary), faces);
            return std::forward<decltype(self)>(self);
        }

        template<boundary::IsBoundary B>
            requires traits::template is_valid_boundary_v<B>
        auto&& with_boundaries(this auto&& self, const std::array<B, 6>& boundaries) {
            self.set_boundaries(boundaries);
            return std::forward<decltype(self)>(self);
        }

        template<controller::IsController C>
            requires traits::template is_valid_controller_v<std::remove_cvref_t<C>>
        auto&& with_controller(this auto&& self, C&& controller) {
            self.add_controller(std::forward<C>(controller));
            return std::forward<decltype(self)>(self);
        }

        template<controller::IsController... C>
        auto&& with_controllers(this auto&& self, C&&... controllers) {
            self.add_controllers(std::forward<C>(controllers)...);
            return std::forward<decltype(self)>(self);
        }

        template<field::IsField F>
            requires traits::template is_valid_field_v<std::remove_cvref_t<F>>
        auto&& with_field(this auto&& self, F&& field) {
            self.add_field(std::forward<F>(field));
            return std::forward<decltype(self)>(self);
        }

        template<field::IsField... F>
        auto&& with_fields(this auto&& self, F&&... fields) {
            self.add_fields(std::forward<F>(fields)...);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_origin(this auto&& self, const vec3& o) {
            self.set_origin(o);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_extent(this auto&& self, const vec3& e) {
            self.set_extent(e);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_origin(this auto&& self, double x, double y, double z) {
            self.set_origin(x, y, z);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_extent(this auto&& self, double x, double y, double z) {
            self.set_extent(x, y, z);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_domain(this auto&& self, const Domain& domain) {
            self.set_domain(domain);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_domain_padding(this auto&& self, double margin) {
            self.domain_padding(margin);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_domain_padding(this auto&& self, const vec3& margin) {
            self.domain_padding(margin);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_domain_padding_factor(this auto&& self, double factor) {
            self.domain_padding_factor(factor);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_domain_padding_factor(this auto&& self, const vec3& factor) {
            self.domain_padding_factor(factor);
            return std::forward<decltype(self)>(self);
        }
    };


    /**
     * @brief CTAD Guide to deduce template parameters from constructor arguments.
     * Allows: auto env = Environment(forces<LJ>, boundaries<Reflective>);
     */
    template<class... Args>
    Environment(Args...)
        -> Environment<
            core::internal::get_pack_t<interactions::internal::ForcePack, Args...>,
            core::internal::get_pack_t<boundary::internal::BoundaryPack, Args...>,
            core::internal::get_pack_t<controller::internal::ControllerPack,Args...>,
            core::internal::get_pack_t<field::internal::FieldPack, Args...>,
            core::internal::get_particle_attributes_t<Args...>
        >;


    namespace core {
        namespace internal {
            template<typename T>
            inline constexpr bool is_environment_v = false;

            template<
                interactions::internal::IsForcePack FPack,
                boundary::internal::IsBoundaryPack BPack,
                controller::internal::IsControllerPack CPack,
                field::internal::IsFieldPack FFPack,
                particle::IsParticleAttributes ParticleData
            >
            inline constexpr bool is_environment_v<Environment<FPack, BPack, CPack, FFPack, ParticleData>> = true;
        }

        template<typename T>
        concept IsEnvironment = internal::is_environment_v<std::remove_cvref_t<T>>;
    }
}



