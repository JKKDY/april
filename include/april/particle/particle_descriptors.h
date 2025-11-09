#pragma once
#include <any>

#include "april/common.h"
#include "april/particle/particle_defs.h"


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
    };



}