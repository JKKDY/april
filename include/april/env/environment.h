#pragma once

#include <vector>
#include <limits>
#include <any>

#include "april/common.h"
#include "april/env/domain.h"
#include "april/forces/force.h"
#include "april/boundaries/boundary.h"
#include "april/controllers/controller.h"
#include "april/fields/field.h"
#include "april/env/traits.h"
#include "april/particle/particle_descriptors.h"
#include "april/particle/particle_defs.h"


namespace april::env {
    inline const auto MARGIN_DONT_CARE = vec3(std::numeric_limits<double>::max());


    struct to_type { ParticleType type; };

    struct between_types { ParticleType t1, t2; };

    struct between_ids { ParticleID id1, id2; };



    template<
    force::IsForcePack FPack,
    boundary::IsBoundaryPack BPack,
    controller::IsControllerPack CPack,
    field::IsFieldPack FFPack,
    IsUserData ParticleData
    >
    class Environment {
    public:
        using traits = internal::EnvironmentTraits<FPack, BPack, CPack, FFPack, ParticleData>;

        explicit Environment(FPack, BPack, CPack, FFPack, ParticleData) {}

        // empty convenience constructor
        Environment()
        : Environment(force::forces<>, boundary::boundaries<>, controller::controllers<>, field::fields<>, NoUserData{}) {}

        // accepts any subset & order of packs
        template<class... Args>
        requires (internal::is_any_pack_v<std::remove_cvref_t<Args>> && ...)
        explicit Environment(Args&&...)
            : Environment(
                internal::get_pack_t<force::ForcePack, Args...>{},
                internal::get_pack_t<boundary::BoundaryPack, Args...>{},
                internal::get_pack_t<controller::ControllerPack,Args...>{},
                internal::get_pack_t<field::FieldPack, Args...>{},
                internal::get_user_data_t<Args...>{}
            ) {}

    private:
        typename traits::environment_data_t data;

        friend auto internal::get_env_data<FPack, BPack, CPack, FFPack> (const Environment& env);

    public:
        // --- Add particles ---

        // Single particle
        void add_particle(const Particle& particle) {
            internal::add_particle_impl(data, particle);
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

        // Cuboid
        std::vector<ParticleID> add_particles(const ParticleCuboid& cuboid) {
            return internal::add_cuboid_particles_impl(data, cuboid);
        }

        // Sphere
        std::vector<ParticleID> add_particles(const ParticleSphere& sphere) {
            return internal::add_sphere_particles_impl(data, sphere);
        }

        // --- Add Forces ---
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

        // --- Add Boundaries ---
        // Single boundary on one face
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundary(B boundary, const boundary::Face face) {
            data.boundaries[face_to_int(face)].template emplace<B>(std::move(boundary));
        }

        // Same boundary applied to multiple faces
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundaries(B boundary, const std::vector<boundary::Face> & faces) {
            for (const boundary::Face face : faces) {
                data.boundaries[face_to_int(face)].template emplace<B>(boundary);
            }
        }

        // Boundaries provided as array (per-face)
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundaries(const std::array<B, 6> & boundaries) {
            for (const boundary::Face face : boundary::all_faces) {
                data.boundaries[face_to_int(face)].template emplace<B>(boundaries[face_to_int(face)]);
            }
        }

        // --- Add Controllers ---
        template<controller::IsController C> requires traits::template is_valid_controller_v<C>
        void add_controller(C controller) {
            data.controllers.add(controller);
        }

        template<controller::IsController... C>
        void add_controllers(C&&... controllers) {
            (data.controllers.add(std::forward<C>(controllers)), ...);
        }

        // --- Add Fields ---
        template<field::IsField F>  requires traits::template is_valid_field_v<F>
        void add_field(F field) {
            data.fields.add(field);
        }

        template<field::IsField... F>
        void add_fields(F&&... fields) {
            (data.fields.add(std::forward<F>(fields)), ...);
        }

        // --- Set Domain ---
        void set_origin(const vec3& origin) {
            this->data.domain.origin = origin;
        }
        void set_origin(const double x, const double y, const double z) {
            set_origin({x,y,z});
        }

        void set_extent(const vec3& extent) {
            this->data.domain.extent = extent;
        }
        void set_extent(const double x, const double y, const double z) {
            set_extent({x,y,z});
        }

        void set_domain(const Domain& domain) {
            data.domain = domain;
        }

        void auto_domain(const vec3& margin_abs) {
            data.margin_abs = margin_abs;
        }
        void auto_domain(const double margin_abs) {
            auto_domain(vec3{margin_abs});
        }

        void auto_domain_factor(const vec3& margin_fac) {
            data.margin_fac = margin_fac;
        }
        void auto_domain_factor(const double margin_fac) {
            auto_domain(vec3{margin_fac});
        }


        // --- DSL-style chaining helpers ---
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

        Environment& with_particles(const ParticleCuboid& cuboid) {
            add_particles(cuboid);
            return *this;
        }

        Environment& with_particles(const ParticleSphere& sphere) {
            add_particles(sphere);
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
        Environment& with_boundary(B&& boundary, boundary::Face face) {
            set_boundary(std::forward<B>(boundary), face);
            return *this;
        }

        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        Environment& with_boundaries(B&& boundary, const std::vector<boundary::Face>& faces) {
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


    // one CTAD guide to deduce the four template parameters from any-order args
    template<class... Args>
    Environment(Args...)
        -> Environment<
            internal::get_pack_t<force::ForcePack, Args...>,
            internal::get_pack_t<boundary::BoundaryPack, Args...>,
            internal::get_pack_t<controller::ControllerPack,Args...>,
            internal::get_pack_t<field::FieldPack, Args...>,
            internal::get_user_data_t<Args...>
        >;



    template<typename T>
    inline constexpr bool is_environment_v = false;

    template<
        force::IsForcePack FPack,
        boundary::IsBoundaryPack BPack,
        controller::IsControllerPack CPack,
        field::IsFieldPack FFPack,
        IsUserData ParticleData
    >
    inline constexpr bool is_environment_v<Environment<FPack, BPack, CPack, FFPack, ParticleData>> = true;

    template<typename T>
    concept IsEnvironment = is_environment_v<std::remove_cvref_t<T>>;
}
