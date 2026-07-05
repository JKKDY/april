/**
* @file environment.hpp
 * @brief User-facing staging area for simulation configuration.
 * * The Environment class collects particles, interactions,
 * boundary conditions, and global controllers. It acts as a high-level
 * blueprint that is later "lowered" into an optimized april::System.
 */

#pragma once

#include <vector>
#include <any>

#include "april/base/types.hpp"
#include "april/core/domain.hpp"
#include "april/core/internal/environment_traits.hpp"

#include "april/interactions/force.hpp"
#include "april/boundaries/boundary.hpp"
#include "april/controllers/controller.hpp"
#include "april/fields/field.hpp"

#include "april/particle/generators.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/attributes.hpp"


namespace april {

    // scoping markers for interactions
    /// @brief Target a single particle type for self-interaction.
    struct to_type { ParticleType type; };
    /// @brief Target interactions between two particle types.
    struct between_types { ParticleType t1, t2; };
    /// @brief Target interactions between two specific particles.
    struct between_ids { ParticleID id1, id2; };

    /**
      * @brief Container for simulation metadata and initial state.
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

        /// @brief Construct an environment with specific type packs.
        explicit Environment(FPack, BPack, CPack, FFPack, ParticleData) {}

        /// @brief Convenience constructor for an empty, default environment.
        Environment() : Environment(forces<>, boundaries<>, controllers<>, fields<>, NoParticleAttributes{}) {}

        /**
         * @brief Variadic constructor. Accepts any subset of packs in any order.
         * Uses internal template-metaprogramming to sort arguments into the correct class slots.
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
         const vec3& position, const vec3& velocity, const double mass, const ParticleType type = 0, const std::any& user_data = {}) {
            Particle p;
            p.type = type,
            p.position =  position,
            p.velocity = velocity,
            p.mass = mass,
            p.state = ParticleState::ALIVE,
            p.user_data = user_data,
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
        void add_force(F force, to_type scope) {
            data.type_interactions.emplace_back(scope.type, scope.type, typename traits::force_variant_t{std::move(force)});
        }

        /// @brief Define interaction between two distinct particle types.
        template<interactions::IsForce F> requires traits::template is_valid_force_v<F>
        void add_force(F force, between_types scope) {
            data.type_interactions.emplace_back(scope.t1, scope.t2, typename traits::force_variant_t{std::move(force)});
        }

        /// @brief Define interaction between two specific particle IDs (Bonds).
        template<interactions::IsForce F> requires traits::template is_valid_force_v<F>
        void add_force(F force, between_ids scope) {
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
        template<controller::IsController C> requires traits::template is_valid_controller_v<C>
        void add_controller(C controller) {
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
        template<field::IsField F>  requires traits::template is_valid_field_v<F>
        void add_field(F field) {
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
        void auto_domain(const vec3d& margin_abs) {
            data.margin_abs = margin_abs;
        }
        void auto_domain(const double margin_abs) {
            data.margin_abs = vec3d{margin_abs};
        }

        /**
         * @brief Define a relative margin based on the current particle distribution size.
         * @details For example, a factor of 0.5 adds 50% of the particle distribution's
         * width as padding to each side.
         * @note If both margin_abs and margin_fac are set, the larger of the two is used per axis.
         */
        void auto_domain_factor(const vec3d& margin_fac) {
            auto_domain(margin_fac);
        }

        void auto_domain_factor(const double margin_fac) {
            auto_domain(vec3d{margin_fac});
        }


        //-------------------
        // DSL-STYLE CHAINING
        //-------------------
        Environment& with_particle(const Particle& p) {
            add_particle(p);
            return *this;
        }

        Environment& with_particle (
        const vec3& position, const vec3& velocity, const double mass, const ParticleType type = 0, const std::optional<ParticleID> id = {}) {
            add_particle(position, velocity, mass, type, id);
            return *this;
        }

        Environment& with_particles(const std::vector<Particle>& ps) {
            add_particles(ps);
            return *this;
        }

        template<IsParticleGenerator G>
        Environment& with_particles(const G& particles) {
            add_particles(particles);
            return *this;
        }

        template<interactions::IsForce F> requires traits::template is_valid_force_v<F>
        Environment& with_force(F&& force, to_type scope) {
            add_force(std::forward<F>(force), scope);
            return *this;
        }

        template<interactions::IsForce F> requires traits::template is_valid_force_v<F>
        Environment& with_force(F&& force, between_types scope) {
            add_force(std::forward<F>(force), scope);
            return *this;
        }

        template<interactions::IsForce F> requires traits::template is_valid_force_v<F>
        Environment& with_force(F&& force, between_ids scope) {
            add_force(std::forward<F>(force), scope);
            return *this;
        }

        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        Environment& with_boundary(B&& boundary, DomainFace face) {
            set_boundary(std::forward<B>(boundary), face);
            return *this;
        }

        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        Environment& with_boundaries(B&& boundary, const std::vector<DomainFace>& faces) {
            set_boundaries(std::forward<B>(boundary), faces);
            return *this;
        }

        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        Environment& with_boundaries(const std::array<B, 6>& boundaries) {
            set_boundaries(boundaries);
            return *this;
        }

        template<controller::IsController C> requires traits::template is_valid_controller_v<C>
        Environment& with_controller(C controller) {
            add_controller(controller);
            return *this;
        }

        template<controller::IsController... C>
        Environment& with_controllers(C&&... controllers) {
            add_controllers(std::forward<C>(controllers)...);
            return *this;
        }

        template<field::IsField F> requires traits::template is_valid_field_v<F>
        Environment& with_field(F field) {
            add_field(field);
            return *this;
        }

        template<field::IsField... F>
        Environment& with_fields(F&&... fields) {
            add_fields(std::forward<F>(fields)...);
            return *this;
        }

        Environment& with_origin(const vec3& o) {
            set_origin(o);
            return *this;
        }

        Environment& with_extent(const vec3& e) {
            set_extent(e);
            return *this;
        }

        Environment& with_origin(const double x,const double y, const double z) {
            set_origin(x,y,z);
            return *this;
        }

        Environment& with_extent(const double x, const double y, const double z) {
            set_extent(x,y,z);
            return *this;
        }

        Environment& with_domain(const Domain& domain) {
            set_domain(domain);
            return *this;
        }

        Environment& with_auto_domain(const double margin) {
            auto_domain(margin);
            return *this;
        }

        Environment& with_auto_domain(const vec3& margin) {
            auto_domain(margin);
            return *this;
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



