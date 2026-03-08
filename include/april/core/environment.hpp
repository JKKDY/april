#pragma once

#include <vector>
#include <any>

#include "april/base/types.hpp"
#include "april/core/domain.hpp"
#include "april/core/internal/environment_traits.hpp"

#include "april/forces/force.hpp"
#include "april/boundaries/boundary.hpp"
#include "april/controllers/controller.hpp"
#include "april/fields/field.hpp"

#include "april/particle/generators.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/attributes.hpp"


namespace april {

    struct to_type { ParticleType type; };
    struct between_types { ParticleType t1, t2; };
    struct between_ids { ParticleID id1, id2; };

    template<
        force::internal::IsForcePack FPack,
        boundary::internal::IsBoundaryPack BPack,
        controller::internal::IsControllerPack CPack,
        field::internal::IsFieldPack FFPack,
        particle::IsParticleAttributes ParticleData>
    class Environment {
    public:
        using traits = core::internal::EnvironmentTraits<FPack, BPack, CPack, FFPack, ParticleData>;

        explicit Environment(FPack, BPack, CPack, FFPack, ParticleData) {}

        // empty convenience constructor
        Environment()
        : Environment(forces<>, boundaries<>, controllers<>, fields<>, NoParticleAttributes{}) {}

        // accepts any subset & order of packs
        template<class... Args>
        requires (core::internal::is_any_pack_v<std::remove_cvref_t<Args>> && ...) &&
            (!std::same_as<std::remove_cvref_t<Args>, Environment> && ...) // make sonarqube shut up about perfect forwarding
        explicit Environment(Args&&...)
            : Environment(
                core::internal::get_pack_t<force::internal::ForcePack, Args...>{},
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

        // Single particle
        void add_particle(const Particle& particle) {
            if (particle.id.has_value() && data.user_particle_ids.contains(particle.id.value())) {
                throw std::invalid_argument("specified id is not unique");
            }

            data.particles.push_back(particle);

            if (particle.id.has_value()) {
                data.user_particle_ids.insert(particle.id.value());
            }

            data.user_particle_types.insert(particle.type);
        }

        // Single particle
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

        // Multiple particles
        void add_particles(const std::vector<Particle>& particles) {
            this->data.particles.reserve(this->data.particles.size() + particles.size());
            for (auto & p : particles)  add_particle(p);
        }

        // Multi particle Initializer
        template<IsParticleGenerator T>
        void add_particles(const T& generator) {
            add_particles(generator.to_particles());
        }


        //-----------------
        // ADD INTERACTIONS
        //-----------------
        // Force applied to a single particle type (self-interaction)
        template<force::IsForce F> requires traits::template is_valid_force_v<F>
        void add_force(F force, to_type scope) {
            // add_force({true, std::pair{scope.type, scope.type}})
            data.type_interactions.emplace_back(scope.type, scope.type, typename traits::force_variant_t{std::move(force)});
        }

        // Force applied between two particle types
        template<force::IsForce F> requires traits::template is_valid_force_v<F>
        void add_force(F force, between_types scope) {
            data.type_interactions.emplace_back(scope.t1, scope.t2, typename traits::force_variant_t{std::move(force)});
        }

        // Force applied between two specific particle IDs
        template<force::IsForce F> requires traits::template is_valid_force_v<F>
        void add_force(F force, between_ids scope) {
            data.id_interactions.emplace_back(scope.id1, scope.id2, typename traits::force_variant_t{std::move(force)});
        }


        //---------------
        // ADD BOUNDARIES
        //---------------
        // Single boundary on one face
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundary(B boundary, const DomainFace face) {
            data.boundaries[boundary::face_to_int(face)].template emplace<B>(std::move(boundary));
        }

        // Same boundary applied to multiple faces
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundaries(B boundary, const std::vector<DomainFace> & faces) {
            for (const DomainFace face : faces) {
                data.boundaries[boundary::face_to_int(face)].template emplace<B>(boundary);
            }
        }

        // Boundaries provided as array (per-face)
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundaries(const std::array<B, 6> & boundaries) {
            for (const DomainFace face : all_faces) {
                data.boundaries[boundary::face_to_int(face)].template emplace<B>(boundaries[boundary::face_to_int(face)]);
            }
        }


        //----------------
        // ADD CONTROLLERS
        //----------------
        template<controller::IsController C> requires traits::template is_valid_controller_v<C>
        void add_controller(C controller) {
            data.controllers.add(controller);
        }

        template<controller::IsController... C>
        void add_controllers(C&&... controllers) {
            (data.controllers.add(std::forward<C>(controllers)), ...);
        }



        //-----------
        // ADD FIELDS
        //-----------
        template<field::IsField F>  requires traits::template is_valid_field_v<F>
        void add_field(F field) {
            data.fields.add(field);
        }

        template<field::IsField... F>
        void add_fields(F&&... fields) {
            (data.fields.add(std::forward<F>(fields)), ...);
        }

        //-----------
        // SET DOMAIN
        //-----------
        void set_origin(const vec3d& origin) {
            this->data.domain.origin = origin;
        }
        void set_origin(const double x, const double y, const double z) {
            set_origin({x,y,z});
        }

        void set_extent(const vec3d& extent) {
            this->data.domain.extent = extent;
        }
        void set_extent(const double x, const double y, const double z) {
            set_extent({x,y,z});
        }

        void set_domain(const Domain& domain) {
            data.domain = domain;
        }

        void auto_domain(const vec3d& margin_abs) {
            data.margin_abs = margin_abs;
        }
        void auto_domain(const double margin_abs) {
            auto_domain(vec3d{margin_abs});
        }

        void auto_domain_factor(const vec3d& margin_fac) {
            data.margin_fac = margin_fac;
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


        template<force::IsForce F> requires traits::template is_valid_force_v<F>
        Environment& with_force(F&& force, to_type scope) {
            add_force(std::forward<F>(force), scope);
            return *this;
        }

        template<force::IsForce F> requires traits::template is_valid_force_v<F>
        Environment& with_force(F&& force, between_types scope) {
            add_force(std::forward<F>(force), scope);
            return *this;
        }

        template<force::IsForce F> requires traits::template is_valid_force_v<F>
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


    // one CTAD guide to deduce the five template parameters from any-order args
    template<class... Args>
    Environment(Args...)
        -> Environment<
            core::internal::get_pack_t<force::internal::ForcePack, Args...>,
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
                force::internal::IsForcePack FPack,
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















