#pragma once

#include <vector>
#include <functional>
#include <unordered_set>
#include <ranges>
#include <algorithm>
#include <variant>
#include <limits>

#include "april/common.h"
#include "april/env/domain.h"
#include "april/env/particle.h"
#include "april/forces/force.h"
#include "april/boundaries/boundary.h"
#include "april/boundaries/absorb.h"

namespace april::env {
    template<class FPack, class BPack> class Environment;

    inline const auto EXTENT_NOT_SET = vec3(std::numeric_limits<double>::max());
    inline const auto ORIGIN_NOT_SET = vec3(std::numeric_limits<double>::max());
    inline const auto MARGIN_DONT_CARE = vec3(std::numeric_limits<double>::max());
    inline const auto ZERO_THERMAL_V = [](const Particle&) {return vec3{}; };


    namespace internal {
        template<force::internal::ForceVariant FV, boundary::internal::BoundaryVariant BV>
        struct EnvironmentData {
            using force_variant_t = FV;
            using boundary_variant_t = BV;

            Domain domain = {ORIGIN_NOT_SET, EXTENT_NOT_SET};
            vec3 margin_abs = {0, 0, 0};
            vec3 margin_fac = {0.5, 0.5, 0.5}; // 50 % margin on each side by default

            std::unordered_set<env::ParticleID> usr_particle_ids;
            std::unordered_set<env::ParticleType> usr_particle_types;

            std::vector<env::Particle> particles;
            std::vector<force::internal::InteractionInfo<force_variant_t>> interactions {};
            std::array<boundary_variant_t, 6> boundaries;
        };

        template<class FPack, class BPack> auto& get_env_data(Environment<FPack, BPack>& env) {
            for (auto & v : env.data.boundaries) {
                if (std::holds_alternative<std::monostate>(v)) {
                    v.template emplace<boundary::Absorb>(); // default-construct Absorb
                }
            }

            return env.data;
        }

    }


    struct ParticleCuboid {
        vec3 origin;
        vec3 mean_velocity;
        uint3 particle_count;
        double distance;
        double particle_mass;
        int type_id;
        std::function<vec3(const Particle&)> thermal_velocity = ZERO_THERMAL_V;
        ParticleState particle_state = ParticleState::ALIVE;

        // fluent setters
        [[nodiscard]] ParticleCuboid& at(const vec3& p) noexcept;
        [[nodiscard]] ParticleCuboid& velocity(const vec3& v) noexcept;
        [[nodiscard]] ParticleCuboid& count(const uint3& n) noexcept;
        [[nodiscard]] ParticleCuboid& spacing(double d) noexcept;
        [[nodiscard]] ParticleCuboid& mass(double m) noexcept;
        [[nodiscard]] ParticleCuboid& type(int t) noexcept;
        [[nodiscard]] ParticleCuboid& thermal(std::function<vec3(const Particle&)> tv);
        [[nodiscard]] ParticleCuboid& state(ParticleState s) noexcept;
    };


    struct ParticleSphere {
        vec3 center;
        vec3 mean_velocity;
        vec3 radii;  // for true sphere set all equal
        double distance;  // packing spacing
        double particle_mass;
        int type_id;
        std::function<vec3(const Particle&)> thermal_velocity = ZERO_THERMAL_V;
        ParticleState particle_state = ParticleState::ALIVE;

        // fluent setters
        [[nodiscard]] ParticleSphere& at(const vec3& c) noexcept;
        [[nodiscard]] ParticleSphere& velocity(const vec3& v) noexcept;
        [[nodiscard]] ParticleSphere& radius_xyz(const vec3& r) noexcept;
        [[nodiscard]] ParticleSphere& radius(double r) noexcept;   // convenience: uniform
        [[nodiscard]] ParticleSphere& spacing(double d) noexcept;
        [[nodiscard]] ParticleSphere& mass(double m) noexcept;
        [[nodiscard]] ParticleSphere& type(int t) noexcept;
        [[nodiscard]] ParticleSphere& thermal(std::function<vec3(const Particle&)> tv);
        [[nodiscard]] ParticleSphere& state(ParticleState s) noexcept;
    };


    struct to_type {
        ParticleType type;
    };

    struct between_types {
        ParticleType t1, t2;
    };

    struct between_ids {
        ParticleID id1, id2;
    };


    template<class... Fs, class... BCs>
    class Environment<force::ForcePack<Fs...>, boundary::BoundaryPack<BCs...>>{
    public:
        explicit Environment(force::ForcePack<Fs...>, boundary::BoundaryPack<BCs...>) {}

        explicit Environment(force::ForcePack<Fs...> force_types)
            : Environment(force_types, boundary::boundaries<>) {}

        static constexpr bool has_absorb = (std::is_same_v<boundary::Absorb, BCs> || ...);

        // expose types for deduction by ForceManager/build_system
        using forces_pack_t = force::ForcePack<Fs...>;
        using boundary_pack_t = boundary::BoundaryPack<BCs...>;

        using force_variant_t = std::variant<Fs...>;
        using boundary_variant_t =
            std::conditional_t<has_absorb,
            std::variant<std::monostate, BCs...>,   // already contains Absorb
            std::variant<std::monostate, boundary::Absorb, BCs...>    // ensure Absorb is present
        >;

    private:
        internal::EnvironmentData<force_variant_t, boundary_variant_t> data;
        template<class FPack, class BPack> friend auto& internal::get_env_data(Environment<FPack, BPack>& env);

    public:

        // --- Add particles ---

        // Single particle
        void add_particle(const Particle& particle) {
            if (particle.id != PARTICLE_ID_DONT_CARE && data.usr_particle_ids.contains(particle.id)) {
                throw std::invalid_argument("specified id is not unique");
            }

            data.particles.push_back(particle);

            if (particle.id != PARTICLE_ID_DONT_CARE) {
                data.usr_particle_ids.insert(particle.id);
            }

            data.usr_particle_types.insert(particle.type);
        }

        // Single particle
        void add_particle(
         const vec3& position, const vec3& velocity, const double mass, const ParticleType type = 0, const ParticleID id = PARTICLE_ID_DONT_CARE) {
            add_particle(Particle{
               .id = id,
               .type = type,
               .position =  position,
               .velocity = velocity,
               .mass = mass,
               .state = ParticleState::ALIVE
           });
        }

        // Multiple particles
        void add_particles(const std::vector<Particle>& particles) {
            this->data.particles.reserve(this->data.particles.size() + particles.size());

            for (auto & p : particles) {
                add_particle(p);
            }
        }

        // Cuboid
        std::vector<ParticleID> add_particles(const ParticleCuboid& cuboid) {
            const uint32_t particle_count = cuboid.particle_count[0] * cuboid.particle_count[1] * cuboid.particle_count[2];
            const double width = cuboid.distance;

            std::vector<ParticleID> ids;
            ids.reserve(particle_count);

            const auto it = std::ranges::max_element(
                data.particles,
                {},               // default `<` comparator
                &Particle::id       // project each Particle to its `id`
            );
            int id = ( it == data.particles.end() ? 0 : it->id+1 );

            data.particles.reserve(data.particles.size() + particle_count);

            for (unsigned int x = 0; x < cuboid.particle_count[0]; ++x) {
                for (unsigned int y = 0; y < cuboid.particle_count[1]; ++y) {
                    for (unsigned int z = 0; z < cuboid.particle_count[2]; ++z) {

                        ids.push_back(id);

                        Particle p = {
                            .id = id++,
                            .type = cuboid.type_id,
                            .position = cuboid.origin + vec3(x * width, y * width, z * width),
                            .velocity = cuboid.mean_velocity,
                            .mass = cuboid.particle_mass,
                            .state = cuboid.particle_state,
                        };
                        p.velocity += cuboid.thermal_velocity(p);

                        add_particle(p);
                    }
                }
            }
            return ids;
        }

        // Sphere
        std::vector<ParticleID> add_particles(const ParticleSphere& sphere) {
            const double width = sphere.distance;
            const vec3 & r = sphere.radii;

            std::vector<ParticleID> ids;

            // get the maximum current id
            const auto it = std::ranges::max_element(
                data.particles,
                {},                 // default `<` comparator
                &Particle::id       // project each Particle to its `id`
            );
            int id = (it == data.particles.end() ? 0 : it->id+1);

            for (int x = -static_cast<int>(sphere.radii.x/width); x < static_cast<int>(sphere.radii.x/width); ++x) {
                for (int y = -static_cast<int>(sphere.radii.y/width); y < static_cast<int>(sphere.radii.y/width); ++y) {
                    for (int z = -static_cast<int>(sphere.radii.z/width); z < static_cast<int>(sphere.radii.z/width); ++z) {

                        vec3 pos = {x * width, y * width, z * width};
                        const vec3 pos_sq = pos * pos;

                        // if not in ellipsoid skip
                        if (pos_sq.x/(r.x*r.x) + pos_sq.y/(r.y*r.y) + pos_sq.z/(r.z*r.z) > 1) continue;

                        ids.push_back(id);
                        Particle p =  {
                            .id = id++,
                            .type = sphere.type_id,
                            .position = sphere.center + pos,
                            .velocity = sphere.mean_velocity,
                            .mass = sphere.particle_mass,
                            .state = sphere.particle_state,
                        };
                        p.velocity += sphere.thermal_velocity(p);

                        add_particle(p);
                    }
                }
            }
            return ids;
        }

        // --- Add Forces ---
        // Force applied to all particles of a given type
        template<force::IsForce F> requires same_as_any<F, Fs...>
        void add_force(F force, to_type scope) {
            data.interactions.emplace_back(true, std::pair{scope.type, scope.type}, force_variant_t{std::move(force)});
        }

        // Force applied between two particle types
        template<force::IsForce F> requires same_as_any<F, Fs...>
        void add_force(F force, between_types scope) {
            data.interactions.emplace_back(true, ParticleTypePair{scope.t1, scope.t2}, force_variant_t{std::move(force)});
        }

        // Force applied between two specific particle IDs
        template<force::IsForce F> requires same_as_any<F, Fs...>
        void add_force(F force, between_ids scope) {
            data.interactions.emplace_back(false, ParticleIDPair{scope.id1, scope.id2}, force_variant_t{std::move(force)});
        }

        // --- Add Boundaries ---
        // Single boundary on one face
        template<boundary::IsBoundary B> requires same_as_any<B, BCs...>
        void set_boundary(B boundary, const boundary::Face face) {
            data.boundaries[face_to_int(face)].template emplace<B>(std::move(boundary));
        }

        // Same boundary applied to multiple faces
        template<boundary::IsBoundary B> requires same_as_any<B, BCs...>
        void set_boundaries(B boundary, const std::vector<boundary::Face> & faces) {
            for (const boundary::Face face : faces) {
                data.boundaries[face_to_int(face)].template emplace<B>(boundary);
            }
        }

        // Boundaries provided as array (per-face)
        template<boundary::IsBoundary B> requires same_as_any<B, BCs...>
        void set_boundaries(const std::array<B, 6> & boundaries) {
            for (const boundary::Face face : boundary::all_faces) {
                data.boundaries[face_to_int(face)].template emplace<B>(boundaries[face_to_int(face)]);
            }
        }

        // --- Set Domain ---
        void set_origin(const vec3& origin) { this->data.domain.origin = origin; }
        void set_origin(const double x, const double y, const double z) { set_origin({x,y,z}); }

        void set_extent(const vec3& extent) { this->data.domain.extent = extent; }
        void set_extent(const double x, const double y, const double z) { set_extent({x,y,z}); }

        void set_domain(const Domain& domain) { data.domain = domain; }

        void auto_domain(const vec3& margin_abs) { data.margin_abs = margin_abs; }
        void auto_domain(const double margin_abs) {auto_domain(vec3{margin_abs}); }

        void auto_domain_factor(const vec3& margin_fac) { data.margin_fac = margin_fac; }
        void auto_domain_factor(const double margin_fac) {auto_domain(vec3{margin_fac}); }


        // --- DSL-style chaining helpers ---
        Environment& with_particle(const Particle& p) {
            add_particle(p);
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

        template<force::IsForce F> requires same_as_any<F, Fs...>
        Environment& with_force(F&& force, to_type scope) {
            add_force(std::forward<F>(force), scope);
            return *this;
        }

        template<force::IsForce F> requires same_as_any<F, Fs...>
        Environment& with_force(F&& force, between_types scope) {
            add_force(std::forward<F>(force), scope);
            return *this;
        }

        template<force::IsForce F> requires same_as_any<F, Fs...>
        Environment& with_force(F&& force, between_ids scope) {
            add_force(std::forward<F>(force), scope);
            return *this;
        }

        template<boundary::IsBoundary B> requires same_as_any<B, BCs...>
        Environment& with_boundary(B&& boundary, boundary::Face face) {
            set_boundary(std::forward<B>(boundary), face);
            return *this;
        }

        template<boundary::IsBoundary B> requires same_as_any<B, BCs...>
        Environment& with_boundaries(B&& boundary, const std::vector<boundary::Face>& faces) {
            set_boundaries(std::forward<B>(boundary), faces);
            return *this;
        }

        template<boundary::IsBoundary B> requires same_as_any<B, BCs...>
        Environment& with_boundaries(const std::array<B, 6>& boundaries) {
            set_boundaries(boundaries);
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

        Environment& with_domain(const vec3& margin) {
            auto_domain(margin);
            return *this;
        }
    };




    // Deduction guide: 2 args
    template<force::IsForce... Fs, boundary::IsBoundary... BCs>
    Environment(force::ForcePack<Fs...>, boundary::BoundaryPack<BCs...>)
        -> Environment<force::ForcePack<Fs...>, boundary::BoundaryPack<BCs...>>;

    // Deduction guide: 1 args
    template<force::IsForce... Fs>
    Environment(force::ForcePack<Fs...>)
    -> Environment<force::ForcePack<Fs...>, boundary::BoundaryPack<>>;

}
