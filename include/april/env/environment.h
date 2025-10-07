#pragma once

#include <vector>
#include <functional>
#include <unordered_set>
#include <ranges>
#include <algorithm>
#include <variant>
#include <limits>

#include "april/common.h"
#include "april/env/particle.h"
#include "april/env/force.h"
#include "april/domain/boundary.h"
#include "april/domain/domain.h"

namespace april::env {
    template<class FPack, class BPack> class Environment;

    inline const auto EXTENT_AUTO = vec3(std::numeric_limits<double>::max());
    inline const auto ORIGIN_AUTO = vec3(std::numeric_limits<double>::max());
    inline const auto ZERO_THERMAL_V = [](const Particle&) {return vec3{}; };


    namespace impl {
        template<ForceVariant FV, IsBoundaryVariant BV> struct EnvironmentData {
            using force_variant_t = FV;
            using boundary_variant_t = BV;

            Domain domain = {EXTENT_AUTO, ORIGIN_AUTO};

            std::unordered_set<env::ParticleID> usr_particle_ids;
            std::unordered_set<env::ParticleType> usr_particle_types;

            std::vector<env::Particle> particles;
            std::vector<InteractionInfo<force_variant_t>> interactions {};
            std::array<boundary_variant_t, 6> boundaries;
        };

        template<class FPack, class BPack> auto& get_env_data(Environment<FPack, BPack>& env) {
            for (auto & v : env.data.boundaries) {
                if (std::holds_alternative<std::monostate>(v)) {
                    v.template emplace<Absorb>(); // default-construct Absorb
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
    class Environment<ForcePack<Fs...>, BoundaryPack<BCs...>>{
    public:
        explicit Environment(ForcePack<Fs...>, BoundaryPack<BCs...>) {}

        explicit Environment(ForcePack<Fs...> force_types)
            : Environment(force_types, boundaries<>) {}

        static constexpr bool has_absorb = (std::is_same_v<Absorb, BCs> || ...);

        // expose types for deduction by ForceManager/build_system
        using forces_pack_t = ForcePack<Fs...>;
        using boundary_pack_t = BoundaryPack<BCs...>;

        using force_variant_t = std::variant<Fs...>;
        using boundary_variant_t =
            std::conditional_t<has_absorb,
            std::variant<std::monostate, BCs...>,   // already contains Absorb
            std::variant<std::monostate, Absorb, BCs...>    // ensure Absorb is present
        >;

    private:
        impl::EnvironmentData<force_variant_t, boundary_variant_t> data;
        template<class FPack, class BPack> friend auto& impl::get_env_data(Environment<FPack, BPack>& env);

    public:
        void add(
            const vec3& position, const vec3& velocity, const double mass, const ParticleType type = 0, const ParticleID id = PARTICLE_ID_DONT_CARE) {
            add(Particle{
               .id = id,
               .type = type,
               .position =  position,
               .velocity = velocity,
               .mass = mass,
               .state = ParticleState::ALIVE
           });
        }

        void add(const Particle& particle) {
            if (particle.id != PARTICLE_ID_DONT_CARE && data.usr_particle_ids.contains(particle.id)) {
                throw std::invalid_argument("specified id is not unique");
            }

            data.particles.push_back(particle);

            if (particle.id != PARTICLE_ID_DONT_CARE) {
                data.usr_particle_ids.insert(particle.id);
            }

            data.usr_particle_types.insert(particle.type);
        }

        void add(const std::vector<Particle>& particles) {
            this->data.particles.reserve(this->data.particles.size() + particles.size());

            for (auto & p : particles) {
                add(p);
            }
        }

        std::vector<ParticleID> add(const ParticleCuboid& cuboid) {
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

                        add(p);
                    }
                }
            }
            return ids;
        }

        std::vector<ParticleID> add(const ParticleSphere& sphere) {
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
                        const vec3 pos_sq = pos.mul(pos);

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

                        add(p);
                    }
                }
            }
            return ids;
        }

        template<IsForce F> requires same_as_any<F, Fs...>
        void add_force(F force, to_type scope) {
            data.interactions.emplace_back(true, std::pair{scope.type, scope.type}, force_variant_t{std::move(force)});
        }

        template<IsForce F> requires same_as_any<F, Fs...>
        void add_force(F force, between_types scope) {
            data.interactions.emplace_back(true, ParticleTypePair{scope.t1, scope.t2}, force_variant_t{std::move(force)});
        }
        template<IsForce F> requires same_as_any<F, Fs...>
        void add_force(F force, between_ids scope) {
            data.interactions.emplace_back(false, ParticleIDPair{scope.id1, scope.id2}, force_variant_t{std::move(force)});
        }

        template<IsBoundary B> requires same_as_any<B, BCs...>
        void add_boundary(B boundary, const Face face) {
            data.boundaries[to_int(face)] = boundary;
        }

        template<IsBoundary B> requires same_as_any<B, BCs...>
        void add_boundary(B boundary, const std::vector<Face> & faces) {
            for (const Face face : faces) {
                data.boundaries[to_int(face)] = boundary;
            }
        }

        void set_origin(const vec3& origin) {
            this->data.domain.origin = origin;
        }

        void set_extent(const vec3& extent) {
            if (extent.x < 0 || extent.y < 0 || extent.z < 0) {
                throw std::invalid_argument("Extent cannot be negative. Got: " + extent.to_string());
            }
            this->data.domain.extent = extent;
        }

        void set_domain(const Domain& domain) {
            data.domain = domain;
        }

        void auto_extent(double) {
            throw std::logic_error("Not implemented yet");
        }

        void add_force_field() {
            throw std::logic_error("Not implemented yet");
        }

        void add_barrier() {
            throw std::logic_error("Not implemented yet");
        }

        void set_boundary_conditions() {
            throw std::logic_error("Not implemented yet");
        }
    };


    // Deduction guide: 2 args
    template<IsForce... Fs, IsBoundary... BCs>
    Environment(ForcePack<Fs...>, BoundaryPack<BCs...>)
        -> Environment<ForcePack<Fs...>, BoundaryPack<BCs...>>;

    // Deduction guide: 1 args
    template<IsForce... Fs>
    Environment(ForcePack<Fs...>)
    -> Environment<ForcePack<Fs...>, BoundaryPack<>>;

}
