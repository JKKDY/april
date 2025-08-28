#pragma once

#include <vector>
#include <functional>
#include <unordered_set>

#include "april/env/particle.h"
#include "april/env/interaction.h"

namespace april::env {
    struct Environment;
    inline const auto EXTENT_AUTO = vec3(std::numeric_limits<double>::infinity());
    inline const auto ORIGIN_AUTO = vec3(std::numeric_limits<double>::infinity());
    inline const auto ZERO_THERMAL_V = [](const Particle&) {return vec3{}; };

    struct Domain {
        vec3 extent = {};
        vec3 origin = {};
    };


    namespace impl {
        struct EnvironmentData {
            Domain domain = {EXTENT_AUTO, ORIGIN_AUTO};

            std::unordered_set<env::ParticleID> usr_particle_ids;
            std::unordered_set<env::ParticleType> usr_particle_types;

            std::vector<env::Particle> particles;
            std::vector<InteractionInfo> interactions;
        };

        EnvironmentData & get_env_data(Environment& env);
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

    struct Environment {
        void add(const vec3& position, const vec3& velocity, double mass, ParticleType type=0, ParticleID id = PARTICLE_ID_DONT_CARE);
        void add(const Particle & particle);
        void add(const std::vector<Particle> & particles);

        std::vector<ParticleID> add(const ParticleCuboid& cuboid);
        std::vector<ParticleID> add(const ParticleSphere& sphere);

        template<IsForce F> void add_force(const F& force, to_type scope);
        template<IsForce F> void add_force(const F& force, between_types scope);
        template<IsForce F> void add_force(const F& force, between_ids scope);

        void set_origin(const vec3& origin);
        void set_extent(const vec3& extent);
        void set_domain(const Domain & domain);
        void auto_extent(double margin);

        void add_force_field();
        void add_barrier();
        void set_boundary_conditions();

    private:
        impl::EnvironmentData data;

        friend impl::EnvironmentData & impl::get_env_data(Environment & env);
    };

    template<IsForce F>
    void Environment::add_force(const F& force, to_type scope) {
        auto ptr = std::make_unique<F>(force);
        data.interactions.emplace_back(true, std::pair{scope.type, scope.type}, std::move(ptr));
    }

    template<IsForce F>
    void Environment::add_force(const F& force, between_types scope) {
        auto ptr = std::make_unique<F>(force);
        data.interactions.emplace_back(true, ParticleTypePair{scope.t1, scope.t2}, std::move(ptr));
    }

    template<IsForce F>
    void Environment::add_force(const F& force, between_ids scope) {
        auto ptr = std::make_unique<F>(force);
        data.interactions.emplace_back(false, ParticleIDPair{scope.id1, scope.id2}, std::move(ptr));
    }
}