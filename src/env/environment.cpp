#include "april/env/environment.h"
#include <algorithm>
#include "april/april.h"


namespace april::env {
    // ---- ParticleCuboid ----
    ParticleCuboid& ParticleCuboid::at(const vec3& p) noexcept {
        origin = p; return *this;
    }
    ParticleCuboid& ParticleCuboid::velocity(const vec3& v) noexcept {
        mean_velocity = v; return *this;
    }
    ParticleCuboid& ParticleCuboid::count(const uint3& n) noexcept {
        particle_count = n; return *this;
    }
    ParticleCuboid& ParticleCuboid::spacing(const double d) noexcept {
        distance = d; return *this;
    }
    ParticleCuboid& ParticleCuboid::mass(const double m) noexcept {
        particle_mass = m; return *this;
    }
    ParticleCuboid& ParticleCuboid::type(const int t) noexcept {
        type_idx = t; return *this;
    }
    ParticleCuboid& ParticleCuboid::thermal(std::function<vec3(const vec3&)> tv) {
        thermal_velocity = std::move(tv); return *this;
    }
    ParticleCuboid& ParticleCuboid::state(const ParticleState s) noexcept {
        particle_state = s; return *this;
    }
    ParticleCuboid& ParticleCuboid::with_data(const std::any & data) noexcept {
        user_data = data; return *this;
    }

    // ---- ParticleSphere ----
    ParticleSphere& ParticleSphere::at(const vec3& c) noexcept {
        center = c; return *this;
    }
    ParticleSphere& ParticleSphere::velocity(const vec3& v) noexcept {
        mean_velocity = v; return *this;
    }
    ParticleSphere& ParticleSphere::radius_xyz(const vec3& r) noexcept {
        radii = r; return *this;
    }
    ParticleSphere& ParticleSphere::radius(double r) noexcept {
        radii = {r, r, r}; return *this;
    }
    ParticleSphere& ParticleSphere::spacing(const double d) noexcept {
        distance = d; return *this;
    }
    ParticleSphere& ParticleSphere::mass(const double m) noexcept {
        particle_mass = m; return *this;
    }
    ParticleSphere& ParticleSphere::type(const int t) noexcept {
        type_idx = t; return *this;
    }
    ParticleSphere& ParticleSphere::thermal(std::function<vec3(const vec3&)> tv) {
        thermal_velocity = std::move(tv); return *this;
    }
    ParticleSphere& ParticleSphere::state(const ParticleState s) noexcept {
        particle_state = s; return *this;
    }
    ParticleSphere& ParticleSphere::with_data(const std::any & data) noexcept {
        user_data = data; return *this;
    }
} // namespace april::env

