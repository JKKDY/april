#pragma once

#include <vector>
#include <functional>
#include <unordered_set>
#include <unordered_map>

#include "april/env/particle.h"
#include "april/env/interaction.h"
#include "april/utils/vec3.hpp"


namespace april::env {
    namespace impl {
        class ParticleIterator;
    }

    enum class Dimension {
        TWO = 2,
        THREE = 3,
        INFER = -1
    };

    struct ParticleCuboid {
        vec3 origin;
        vec3 mean_velocity;
        uint3 particle_count;
        double distance;
        double mass;
        int type;
        std::function<vec3(const Particle&)> thermal_velocity = [](const Particle&) {return vec3{}; };
        ParticleState state = ParticleState::ALIVE;
    };

    struct ParticleSphere {
        vec3 center;
        vec3 mean_velocity;
        vec3 radius;
        double distance;
        double mass;
        int type;
        std::function<vec3(const Particle&)> thermal_velocity = [](const Particle&) {return vec3{}; };
        ParticleState state = ParticleState::ALIVE;
    };


    class Environment {
    public:
        Environment();

        void add_particle(const Particle & particle);
        void add_particles(const std::vector<Particle> & particles);

        std::vector<ParticleID> add_particle_cuboid(const ParticleCuboid& cuboid);
        std::vector<ParticleID> add_particle_sphere(const ParticleSphere& sphere);

        template<IsForce F> void add_force(const F & force, ParticleType type);
        template<IsForce F> void add_force(const F & force, ParticleTypePair types);
        template<IsForce F> void add_interaction(const F & force, ParticleIDPair ids);

        void build();
  
        impl::ParticleIterator particles(ParticleState state = ParticleState::ALL);

    private:
        void validate_parameters();
        void map_ids_and_types_to_internal();
        void build_particles();

        std::vector<Particle> particle_infos;

        std::vector<impl::Particle> particle_storage;
        std::vector<impl::InteractionInfo> interactions;
        impl::InteractionManager interaction_manager;

        bool is_built;
        std::unordered_set<ParticleType> usr_particle_types;
        std::unordered_map<ParticleType, impl::ParticleType> usr_types_to_impl_types;
        std::unordered_set<ParticleID> usr_particle_ids; 
        std::unordered_map<ParticleID, impl::ParticleID> usr_ids_to_impl_ids;
    };


    template<IsForce F> void Environment::add_force(const F & force, ParticleType type) {
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        interactions.emplace_back(true, std::pair{ type, type }, std::move(ptr));
    }
    template<IsForce F> void Environment::add_force(const F & force, ParticleTypePair types) {
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        interactions.emplace_back(true, types, std::move(ptr));
    }
    template<IsForce F> void Environment::add_interaction(const F & force, ParticleIDPair ids) {
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        interactions.emplace_back(false, ids, std::move(ptr));
    }


    namespace impl {

        class ParticleIterator {
            using reference = impl::Particle&;
            using pointer = impl::Particle*;

            class Iterator {
            public:
                Iterator(std::vector<impl::Particle> & particles, ParticleState state, size_t idx) : 
                    particle_storage(particles), state(state), idx(idx) {}

                reference operator*() const {
                    return particle_storage[idx];
                } 

                Iterator& operator++() {  
                    if (++idx == particle_storage.size()) return *this;
                    Particle & p = particle_storage[idx];
                    if (static_cast<unsigned int>(p.state) & static_cast<unsigned int>(state)) 
                        return *this;
                    else 
                        return ++*this;
                }

                Iterator operator++(int) { 
                    auto tmp = *this; 
                    ++*this; 
                    return tmp; 
                }

                bool operator==(const Iterator& other)const { 
                    return idx == other.idx; 
                }

            private:
                std::vector<impl::Particle> & particle_storage;
                ParticleState state;
                size_t idx;
            };
        public:
            ParticleIterator(std::vector<impl::Particle> & particles, ParticleState state) : 
                particle_storage(particles), state(state) {}

            Iterator begin() {return Iterator(particle_storage, state, 0);}
            Iterator end() {return Iterator(particle_storage, state, particle_storage.size());}

        private:
            std::vector<impl::Particle> & particle_storage;
            ParticleState state;
        };

    }
}