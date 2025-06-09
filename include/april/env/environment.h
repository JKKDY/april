#pragma once

#include <vector>
#include <functional>
#include <unordered_set>
#include <unordered_map>

#include "april/env/particle.h"
#include "april/env/interaction.h"
#include "april/containers/container.h"

namespace april::env {
    namespace impl {
        class ParticleIterator;
    }

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
        static const vec3 EXTENT_AUTO;
        static const vec3 ORIGIN_AUTO;

        Environment();

        void add_particle(const vec3& position, const vec3& velocity, double mass, ParticleType type=0, ParticleID id = PARTICLE_ID_UNDEFINED);
        void add_particle(const Particle & particle);
        void add_particles(const std::vector<Particle> & particles);

        std::vector<ParticleID> add_particle_cuboid(const ParticleCuboid& cuboid);
        std::vector<ParticleID> add_particle_sphere(const ParticleSphere& sphere);

        template<IsForce F> void add_force_to_type(const F & force, ParticleType type);
        template<IsForce F> void add_force_between_types(const F & force, ParticleType t1, ParticleType t2);
        template<IsForce F> void add_force_between_ids(const F & force, ParticleID id1, ParticleID id2);

        void set_extent(const vec3 & size);
        void set_origin(const vec3 & origin_vec);
        void set_container(std::unique_ptr<core::Container> container_ptr);

        void build();

        void update_forces() const;

        impl::ParticleIterator particles(ParticleState state = ParticleState::ALL);
        const std::vector<impl::Particle> & export_particles();

    private:
        void validate_inputs();
        void map_ids_and_types_to_internal();
        void build_particles();
        void finalize_environment_size();

        std::vector<Particle> particle_infos;

        bool is_built;

        vec3 extent = EXTENT_AUTO;
        vec3 origin = ORIGIN_AUTO;

        vec3 particle_bbox_extent;
        vec3 particle_bbox_origin;

        std::unique_ptr<core::Container> container;

        std::vector<impl::Particle> particle_storage;
        impl::InteractionManager interaction_manager;

        std::vector<impl::InteractionInfo> interactions;
        std::unordered_set<ParticleType> usr_particle_types;
        std::unordered_map<ParticleType, impl::ParticleType> usr_types_to_impl_types;
        std::unordered_set<ParticleID> usr_particle_ids; 
        std::unordered_map<ParticleID, impl::ParticleID> usr_ids_to_impl_ids;
    };


    template<IsForce F> void Environment::add_force_to_type(const F & force, ParticleType type) {
        if (is_built)
            throw std::logic_error("cannot add force. environment has already been built.");
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        interactions.emplace_back(true, std::pair{ type, type }, std::move(ptr));
    }
    template<IsForce F> void Environment::add_force_between_types(const F & force, ParticleType t1, ParticleType t2) {
        if (is_built)
            throw std::logic_error("cannot add force. environment has already been built.");
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        interactions.emplace_back(true, ParticleTypePair{t1, t2}, std::move(ptr));
    }
    template<IsForce F> void Environment::add_force_between_ids(const F & force, ParticleID id1, ParticleID id2) {
        if (is_built)
            throw std::logic_error("cannot add force. environment has already been built.");
        std::unique_ptr<Force> ptr = std::make_unique<F>(force);
        interactions.emplace_back(false, ParticleIDPair{id1, id2}, std::move(ptr));
    }


    namespace impl {

        class ParticleIterator {
            class Iterator {
            public:
                using iterator_category = std::input_iterator_tag;
                using value_type = impl::Particle;
                using difference_type = std::ptrdiff_t;
                using pointer = impl::Particle*;
                using reference = impl::Particle&;

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
                        return ++*this; //TODO use loop instead of recursion
                }

                Iterator operator++(int) {
                    const auto tmp = *this;
                    ++*this; 
                    return tmp; 
                }

                bool operator==(const Iterator& other)const { 
                    return idx == other.idx; 
                }

                bool operator!=(const Iterator& other) const {
                    return !(*this == other);
                }

            private:
                std::vector<impl::Particle> & particle_storage;
                ParticleState state;
                size_t idx;
            };
        public:
            ParticleIterator(std::vector<impl::Particle> & particles, const ParticleState state) :
                particle_storage(particles), state(state) {}

            [[nodiscard]] auto begin() const {return Iterator(particle_storage, state, 0);}
            [[nodiscard]] auto end() const {return Iterator(particle_storage, state, particle_storage.size());}

        private:
            std::vector<impl::Particle> & particle_storage;
            ParticleState state;
        };

    }
}