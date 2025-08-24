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
    inline const auto zero_velocity = [](const Particle&) {return vec3{}; };

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
        double mass;
        int type;
        std::function<vec3(const Particle&)> thermal_velocity = zero_velocity;
        ParticleState state = ParticleState::ALIVE;
    };

    struct ParticleSphere {
        vec3 center;
        vec3 mean_velocity;
        vec3 radius;
        double distance;
        double mass;
        int type;
        std::function<vec3(const Particle&)> thermal_velocity = zero_velocity;
        ParticleState state = ParticleState::ALIVE;
    };

    struct Environment {
        void add_particle(const vec3& position, const vec3& velocity, double mass, ParticleType type=0, ParticleID id = PARTICLE_ID_DONT_CARE);
        void add_particle(const Particle & particle);
        void add_particles(const std::vector<Particle> & particles);

        std::vector<ParticleID> add_particle_cuboid(const ParticleCuboid& cuboid);
        std::vector<ParticleID> add_particle_sphere(const ParticleSphere& sphere);

        template<IsForce F> void add_force_to_type(const F & force, ParticleType type);
        template<IsForce F> void add_force_between_types(const F & force, ParticleType t1, ParticleType t2);
        template<IsForce F> void add_force_between_ids(const F & force, ParticleID id1, ParticleID id2);

        void set_origin(const vec3& origin);
        void set_extent(const vec3& extent);

        void add_force_field();
        void add_barrier();
        void set_boundary_conditions();

    private:
        impl::EnvironmentData data;

        friend impl::EnvironmentData & impl::get_env_data(Environment & env);
    };


    template<IsForce F> void Environment::add_force_to_type(const F & force, ParticleType type) {
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        data.interactions.emplace_back(true, std::pair{ type, type }, std::move(ptr));
    }
    template<IsForce F> void Environment::add_force_between_types(const F & force, ParticleType t1, ParticleType t2) {
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        data.interactions.emplace_back(true, ParticleTypePair{t1, t2}, std::move(ptr));
    }
    template<IsForce F> void Environment::add_force_between_ids(const F & force, ParticleID id1, ParticleID id2) {
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        data.interactions.emplace_back(false, ParticleIDPair{id1, id2}, std::move(ptr));
    }








    // namespace impl {
    //
    //     class ParticleIterator {
    //         class Iterator {
    //         public:
    //             using iterator_category = std::input_iterator_tag;
    //             using value_type = impl::Particle;
    //             using difference_type = std::ptrdiff_t;
    //             using pointer = impl::Particle*;
    //             using reference = impl::Particle&;
    //
    //             Iterator(std::vector<impl::Particle> & particles, const ParticleState state, const size_t index) :
    //                 particle_storage(particles), state(state), idx(index) {
    //
    //                 while (idx < particle_storage.size() && not
    //                     (static_cast<unsigned int>(particle_storage[idx].state) & static_cast<unsigned int>(state))) {
    //                     idx++;
    //                 }
    //             }
    //
    //             reference operator*() const {
    //                 return particle_storage[idx];
    //             }
    //
    //             Iterator& operator++() {
    //                 if (++idx == particle_storage.size()) return *this;
    //                 const Particle & p = particle_storage[idx];
    //                 if (static_cast<bool>(p.state & state))
    //                     return *this;
    //                 else
    //                     return ++*this; //TODO use loop instead of recursion
    //             }
    //
    //             Iterator operator++(int) {
    //                 const auto tmp = *this;
    //                 ++*this;
    //                 return tmp;
    //             }
    //
    //             bool operator==(const Iterator& other)const {
    //                 return idx == other.idx;
    //             }
    //
    //             bool operator!=(const Iterator& other) const {
    //                 return !(*this == other);
    //             }
    //
    //         private:
    //             std::vector<impl::Particle> & particle_storage;
    //             ParticleState state;
    //             size_t idx;
    //         };
    //     public:
    //         ParticleIterator(std::vector<impl::Particle> & particles, const ParticleState state) :
    //             particle_storage(particles), state(state) {}
    //
    //         [[nodiscard]] auto begin() const {return Iterator(particle_storage, state, 0);}
    //         [[nodiscard]] auto end() const {return Iterator(particle_storage, state, particle_storage.size());}
    //
    //     private:
    //         std::vector<impl::Particle> & particle_storage;
    //         ParticleState state;
    //     };
    // }
}