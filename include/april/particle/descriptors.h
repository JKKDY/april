#pragma once
#include <any>

#include "april/common.h"
#include "april/particle/defs.h"


namespace april::env {

    // user facing declaration with optional fields and non typed field for user data
    struct Particle {
        std::optional<ParticleID> id;			// The id of the particle.
        ParticleType type = 0;  				// The type of the particle.

        vec3 position;      					// The position of the particle.
        vec3 velocity;      					// The velocity of the particle.

        double mass{};        					// The mass of the particle.
        ParticleState state{};					// The state of the particle.

        // optional data e.g. if initializing from a simulation snapshot
        std::optional<vec3> old_position;		// previous position of the particle. Useful for applying boundary conditions
        std::optional<vec3> old_force;			// previous force acting on the particle.
        std::optional<vec3> force;				// current force

        std::any user_data {}; // custom user data

        [[nodiscard]] Particle& with_id(ParticleID v) noexcept {
            id = v; return *this;
        }
        [[nodiscard]] Particle& as_type(const ParticleType v) noexcept {
            type = v; return *this;
        }
        [[nodiscard]] Particle& at(const vec3& v) noexcept {
            position = v; return *this;
        }
        [[nodiscard]] Particle& at(const double x, const double y, const double z) noexcept {
            position = {x,y,z}; return *this;
        }
        [[nodiscard]] Particle& with_velocity(const vec3& v) noexcept {
            velocity = v; return *this;
        }
        [[nodiscard]] Particle& with_velocity(const double x, const double y, const double z) noexcept {
            velocity = {x,y,z}; return *this;
        }
        [[nodiscard]] Particle& with_mass(const double v) noexcept {
            mass = v; return *this;
        }
        [[nodiscard]] Particle& with_state(const ParticleState v) noexcept {
            state = v; return *this;
        }
        [[nodiscard]] Particle& with_old_position(const vec3& v) noexcept {
            old_position = v; return *this;
        }
        [[nodiscard]] Particle& with_old_force(const vec3& v) noexcept {
            old_force = v; return *this;
        }
        [[nodiscard]] Particle& with_force(const vec3& v) noexcept {
            force = v; return *this;
        }
        [[nodiscard]] Particle& with_data(const std::any& v) noexcept {
            user_data = v; return *this; }
    };


	inline const auto ZERO_THERMAL_V = [](const vec3&) {return vec3{}; };

    template<typename T>
    concept IsParticleGenerator = requires (T x) {
        {x.to_particles()} -> std::same_as<std::vector<Particle>>;
    };

    struct ParticleCuboid {
        vec3 origin;
        vec3 mean_velocity;
        uint3 particle_count;
        double distance;
        double particle_mass;
        ParticleType type_idx;
        std::any user_data;
        std::function<vec3(const vec3&)> thermal_velocity = ZERO_THERMAL_V;
        ParticleState particle_state = ParticleState::ALIVE;

        // fluent setters
        [[nodiscard]] ParticleCuboid& at(const vec3& p) noexcept {
            origin = p; return *this;
        }
        [[nodiscard]] ParticleCuboid& at(const double x, const double y, const double z) noexcept {
            origin = {x,y,z}; return *this;
        }
        [[nodiscard]] ParticleCuboid& velocity(const vec3& v) noexcept {
            mean_velocity = v; return *this;
        }
        [[nodiscard]] ParticleCuboid& velocity(const double x, const double y, const double z) noexcept {
            mean_velocity = {x,y,z}; return *this;
        }
        [[nodiscard]] ParticleCuboid& count(const uint3& n) noexcept {
            particle_count = n; return *this;
        }
        [[nodiscard]] ParticleCuboid& count(const unsigned x, const unsigned y, const unsigned z) noexcept {
            particle_count = {x,y,z}; return *this;
        }
        [[nodiscard]] ParticleCuboid& spacing(const double d) noexcept {
            distance = d; return *this;
        }
        [[nodiscard]] ParticleCuboid& mass(const double m) noexcept {
            particle_mass = m; return *this;
        }
        [[nodiscard]] ParticleCuboid& type(const int t) noexcept {
            type_idx = t; return *this;
        }
        [[nodiscard]] ParticleCuboid& thermal(std::function<vec3(const vec3&)> tv) {
            thermal_velocity = std::move(tv); return *this;
        }
        [[nodiscard]] ParticleCuboid& state(const ParticleState s) noexcept {
            particle_state = s; return *this;
        }
        [[nodiscard]] ParticleCuboid& with_data(const std::any & data) noexcept {
            user_data = data; return *this;
        }

        [[nodiscard]] std::vector<Particle> to_particles() const {
            if (distance == 0) {
                throw std::logic_error("Cuboid inter-particle distance is set to 0!");
            }
            std::vector<Particle> particles;
            particles.reserve(particle_count.x * particle_count.y * particle_count.z);

            for (unsigned int x = 0; x < particle_count.x; ++x) {
                for (unsigned int y = 0; y < particle_count.y; ++y) {
                    for (unsigned int z = 0; z < particle_count.z; ++z) {
                        Particle p;

                        p.id = std::nullopt;
                        p.type		= type_idx;
                        p.position	= origin + vec3(x * distance, y * distance, z * distance);
                        p.velocity	= mean_velocity;
                        p.mass		= particle_mass;
                        p.state		= particle_state;
                        p.user_data	= user_data;
                        p.velocity += thermal_velocity(p.position);

                        particles.push_back(p);
                    }
                }
            }

            return particles;
        }
    };


    struct ParticleSphere {
        vec3 center;
        vec3 mean_velocity;
        vec3 radii;  // for true sphere set all equal
        double distance;  // packing spacing
        double particle_mass;
        ParticleType type_idx;
        std::any user_data;
        std::function<vec3(const vec3&)> thermal_velocity = ZERO_THERMAL_V;
        ParticleState particle_state = ParticleState::ALIVE;

        // fluent setters
        [[nodiscard]] ParticleSphere& at(const vec3& c) noexcept {
            center = c; return *this;
        }
        [[nodiscard]] ParticleSphere& at(const double x, const double y, const double z) noexcept {
            center = {x,y,z}; return *this;
        }
        [[nodiscard]] ParticleSphere& velocity(const vec3& v) noexcept {
            mean_velocity = v; return *this;
        }
        [[nodiscard]] ParticleSphere& velocity(const double x, const double y, const double z) noexcept {
            mean_velocity = {x,y,z}; return *this;
        }
        [[nodiscard]] ParticleSphere& radius_xyz(const vec3& r) noexcept {
            radii = r; return *this;
        }
        [[nodiscard]] ParticleSphere& radius_xyz(const double x, const double y, const double z) noexcept {
            radii = {x,y,z}; return *this;
        }
        [[nodiscard]] ParticleSphere& radius(double r) noexcept {   // convenience: uniform
            radii = {r, r, r}; return *this;
        }
        [[nodiscard]] ParticleSphere& spacing(const double d) noexcept {
            distance = d; return *this;
        }
        [[nodiscard]] ParticleSphere& mass(const double m) noexcept {
            particle_mass = m; return *this;
        }
        [[nodiscard]] ParticleSphere& type(const int t) noexcept {
            type_idx = t; return *this;
        }
        [[nodiscard]] ParticleSphere& thermal(std::function<vec3(const vec3&)> tv) {
            thermal_velocity = std::move(tv); return *this;
        }
        [[nodiscard]] ParticleSphere& state(const ParticleState s) noexcept {
            particle_state = s; return *this;
        }
        [[nodiscard]] ParticleSphere& with_data(const std::any & data) noexcept {
            user_data = data; return *this;
        }
        [[nodiscard]] std::vector<Particle> to_particles() const {
            if (distance == 0) {
                throw std::logic_error("Sphere inter-particle distance is set to 0");
            }


            const vec3 radii = {
                std::max(radii.x, distance),
                std::max(radii.y, distance),
                std::max(radii.z, distance)
            };
            const double ellipsoid_volume = 4.0/3.0*3.14*radii.x*radii.y*radii.z;

            std::vector<Particle> particles;
            particles.reserve(static_cast<size_t>(ellipsoid_volume / (distance*distance*distance)));

            for (int x = -static_cast<int>(radii.x/distance); x < static_cast<int>(radii.x/distance); ++x) {
                for (int y = -static_cast<int>(radii.y/distance); y < static_cast<int>(radii.y/distance); ++y) {
                    for (int z = -static_cast<int>(radii.z/distance); z < static_cast<int>(radii.z/distance); ++z) {

                        const vec3 pos = {x * distance, y * distance, z * distance};
                        const vec3 pos_sq = pos * pos;

                        // if not in ellipsoid skip
                        if (pos_sq.x/(radii.x*radii.x) +
                            pos_sq.y/(radii.y*radii.y) +
                            pos_sq.z/(radii.z*radii.z) >= 1) continue;

                        Particle p;
                        p.id        = std::nullopt;
                        p.type      = type_idx;
                        p.position  = center + pos;
                        p.velocity  = mean_velocity;
                        p.mass      = particle_mass;
                        p.state     = particle_state;
                        p.user_data = user_data;
                        p.velocity += thermal_velocity(p.position);

                        particles.push_back(p);
                    }
                }
            }

            return particles;
        }
    };



}