#pragma once
#include <any>

#include "april/base/types.hpp"
#include "april/particle/defs.hpp"
#include "april/particle/particle.hpp"

namespace april::env {

    inline const auto ZERO_THERMAL_V = [](const vec3&) {return vec3{}; };

    template<typename T>
    concept IsParticleGenerator = requires (T x) {
        {x.to_particles()} -> std::same_as<std::vector<Particle>>;
    };

    struct ParticleCuboid {
        vec3d origin;
        vec3 mean_velocity;
        uint3 particle_count;
        double distance;
        double particle_mass;
        ParticleType type_idx;
        std::any user_data;
        std::function<vec3(const vec3&)> thermal_velocity = ZERO_THERMAL_V;
        ParticleState particle_state = ParticleState::ALIVE;

        // fluent setters
        ParticleCuboid& at(const vec3& p) noexcept {
            origin = p; return *this;
        }
        ParticleCuboid& at(const vec3::type x, const vec3::type y, const vec3::type z) noexcept {
            origin = {x,y,z}; return *this;
        }
        ParticleCuboid& velocity(const vec3& v) noexcept {
            mean_velocity = v; return *this;
        }
        ParticleCuboid& velocity(const vec3::type x, const vec3::type y, const vec3::type z) noexcept {
            mean_velocity = {x,y,z}; return *this;
        }
        ParticleCuboid& count(const uint3& n) noexcept {
            particle_count = n; return *this;
        }
        ParticleCuboid& count(const unsigned x, const unsigned y, const unsigned z) noexcept {
            particle_count = {x,y,z}; return *this;
        }
        ParticleCuboid& spacing(const double d) noexcept {
            distance = d; return *this;
        }
        ParticleCuboid& mass(const double m) noexcept {
            particle_mass = m; return *this;
        }
        ParticleCuboid& type(const int t) noexcept {
            type_idx = t; return *this;
        }
        ParticleCuboid& thermal(std::function<vec3(const vec3&)> tv) {
            thermal_velocity = std::move(tv); return *this;
        }
        ParticleCuboid& state(const ParticleState s) noexcept {
            particle_state = s; return *this;
        }
        ParticleCuboid& with_data(const std::any & data) noexcept {
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
                        p.type     = type_idx;
                        p.position  = origin + vec3(x * distance, y * distance, z * distance);
                        p.velocity  = mean_velocity;
                        p.mass     = particle_mass;
                        p.state    = particle_state;
                        p.user_data = user_data;
                        p.velocity += thermal_velocity(p.position);

                        particles.push_back(p);
                    }
                }
            }

            return particles;
        }
    };


    struct ParticleSphere {
        vec3 mean_velocity;
        vec3d center;
        vec3d radii;  // for true sphere set all equal
        double distance;  // packing spacing
        double particle_mass;
        ParticleType type_idx;
        std::any user_data;
        std::function<vec3(const vec3&)> thermal_velocity = ZERO_THERMAL_V;
        ParticleState particle_state = ParticleState::ALIVE;

        // fluent setters
        ParticleSphere& at(const vec3& c) noexcept {
            center = c; return *this;
        }
        ParticleSphere& at(const vec3::type x, const vec3::type y, const vec3::type z) noexcept {
            center = {x,y,z}; return *this;
        }
        ParticleSphere& velocity(const vec3& v) noexcept {
            mean_velocity = v; return *this;
        }
        ParticleSphere& velocity(const vec3::type x, const vec3::type y, const vec3::type z) noexcept {
            mean_velocity = {x,y,z}; return *this;
        }
        ParticleSphere& radius_xyz(const vec3& r) noexcept {
            radii = r; return *this;
        }
        ParticleSphere& radius_xyz(const vec3::type x, const vec3::type y, const vec3::type z) noexcept {
            radii = {x,y,z}; return *this;
        }
        ParticleSphere& radius(vec3::type r) noexcept {   // convenience: uniform
            radii = {r, r, r}; return *this;
        }
        ParticleSphere& spacing(const double d) noexcept {
            distance = d; return *this;
        }
        ParticleSphere& mass(const double m) noexcept {
            particle_mass = m; return *this;
        }
        ParticleSphere& type(const int t) noexcept {
            type_idx = t; return *this;
        }
        ParticleSphere& thermal(std::function<vec3(const vec3&)> tv) {
            thermal_velocity = std::move(tv); return *this;
        }
        ParticleSphere& state(const ParticleState s) noexcept {
            particle_state = s; return *this;
        }
        ParticleSphere& with_data(const std::any & data) noexcept {
            user_data = data; return *this;
        }
        [[nodiscard]] std::vector<Particle> to_particles() const {
            if (distance == 0) {
                throw std::logic_error("Sphere inter-particle distance is set to 0");
            }

            const vec3d eff_radii = {
                std::max(radii.x, distance),
                std::max(radii.y, distance),
                std::max(radii.z, distance)
            };
            const double ellipsoid_volume = 4.0/3.0*3.14*eff_radii.x*eff_radii.y*eff_radii.z;

            std::vector<Particle> particles;
            particles.reserve(static_cast<size_t>(ellipsoid_volume / (distance*distance*distance)));

            for (int x = -static_cast<int>(eff_radii.x/distance); x < static_cast<int>(eff_radii.x/distance); ++x) {
                for (int y = -static_cast<int>(eff_radii.y/distance); y < static_cast<int>(eff_radii.y/distance); ++y) {
                    for (int z = -static_cast<int>(eff_radii.z/distance); z < static_cast<int>(eff_radii.z/distance); ++z) {

                        const vec3 pos = vec3{
                            static_cast<vec3::type>(x * distance),
                            static_cast<vec3::type>(y * distance),
                            static_cast<vec3::type>(z * distance)
                        };
                        const vec3 pos_sq = pos * pos;

                        // if not in ellipsoid skip
                        if (pos_sq.x/(eff_radii.x*eff_radii.x) +
                            pos_sq.y/(eff_radii.y*eff_radii.y) +
                            pos_sq.z/(eff_radii.z*eff_radii.z) >= 1) continue;

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