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
        template<IsForce F> void add_force(const F & force, ParticleType t1, ParticleType t2);
        template<IsForce F> void add_interaction(const F & force, ParticleID id1, ParticleID id2);

        void build();
  
        impl::ParticleIterator particles(ParticleState state = ParticleState::ALL);

        const std::vector<impl::Particle> & export_particles();

    private:
        void validate_inputs();
        void map_ids_and_types_to_internal();
        void build_particles();

        void add_interaction(const impl::InteractionInfo & interaction);


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
    template<IsForce F> void Environment::add_force(const F & force, ParticleType t1, ParticleType t2) {
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        interactions.emplace_back(true, ParticleTypePair{t1, t2}, std::move(ptr));
    }
    template<IsForce F> void Environment::add_interaction(const F & force, ParticleID id1, ParticleID id2) {
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        interactions.emplace_back(false, ParticleIDPair{id1, id2}, std::move(ptr));
    }


    namespace impl {

        class ParticleIterator {
            using reference = impl::Particle&;
            using pointer = impl::Particle*;

            class Iterator {
            public:
                Iterator(std::vector<impl::Particle> & particles, const ParticleState state, const size_t index) :
                    particle_storage(particles), state(state), idx(index) {

                    while (idx < particle_storage.size() && not
                        (static_cast<unsigned int>(particle_storage[idx].state) & static_cast<unsigned int>(state))) {
                        idx++;
                    }
                }

                reference operator*() const {
                    return particle_storage[idx];
                } 

                Iterator& operator++() {  
                    if (++idx == particle_storage.size()) return *this;
                    const Particle & p = particle_storage[idx];
                    if (static_cast<bool>(p.state & state))
                        return *this;
                    else 
                        return ++*this;
                }

                Iterator operator++(int) {
                    const auto tmp = *this;
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
            ParticleIterator(std::vector<impl::Particle> & particles, const ParticleState state) :
                particle_storage(particles), state(state) {}

            Iterator begin() const {return Iterator(particle_storage, state, 0);}
            Iterator end() const {return Iterator(particle_storage, state, particle_storage.size());}

        private:
            std::vector<impl::Particle> & particle_storage;
            ParticleState state;
        };

    }
}