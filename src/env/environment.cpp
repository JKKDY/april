#include "april/env/environment.h"
#include <algorithm>
#include <stdexcept>
#include <ranges>

namespace april::env {
    Environment::Environment(): is_built(false) {

    }

    void Environment::add_particle(const Particle &particle) {
        if (particle.id != PARTICLE_ID_UNDEFINED && usr_particle_ids.contains(particle.id)) {
            throw std::invalid_argument("specified id is no unique");
        }

        particle_infos.push_back(particle);   

        // collect id and type information
        if (particle.id != PARTICLE_ID_UNDEFINED) usr_particle_ids.insert(particle.id);
        usr_particle_types.insert(particle.type);
    }


    void Environment::add_particles(const std::vector<Particle> & particles) {
        for (auto & x: particles) {
            add_particle(x);
        }
    }


    std::vector<ParticleID> Environment::add_particle_cuboid(const ParticleCuboid& cuboid) {
        uint32_t particle_count = cuboid.particle_count[0] * cuboid.particle_count[1] * cuboid.particle_count[2];
        double width = cuboid.distance;

        std::vector<ParticleID> ids;
        ids.reserve(particle_count);

        auto it = std::ranges::max_element(
            particle_infos,
            {},                 // default `<` comparator
            &Particle::id       // project each Particle to its `id`
        );
        int id = it == particle_infos.end() ? 0 : it->id+1;

        particle_infos.reserve(particle_infos.size() + particle_count);
        
       for (unsigned int x = 0; x < cuboid.particle_count[0]; ++x) {
            for (unsigned int y = 0; y < cuboid.particle_count[1]; ++y) {
                for (unsigned int z = 0; z < cuboid.particle_count[2]; ++z) {

                    ids.push_back(id);

                    Particle p = {
                        .id = id++,
                        .type = cuboid.type,
                        .position = cuboid.origin + vec3(x * width, y * width, z * width),
                        .velocity = cuboid.mean_velocity,
                        .mass = cuboid.mass,
                        .state = cuboid.state,
                    };
                    p.velocity += cuboid.thermal_velocity(p);

                    add_particle(p);
                }
            }
        }

        return ids;
    }


    std::vector<ParticleID> Environment::add_particle_sphere(const ParticleSphere& sphere) {
        double width = sphere.distance;
        const vec3 & r = sphere.radius; 

        std::vector<ParticleID> ids;

        // get the the maximum current id
        // TODO use projector of projector implciity relies on  PARTICLE_ID_UNDEFINED = min<ParticleID>. Make it explicity e.g. via lambda
        auto it = std::ranges::max_element(
            particle_infos,
            {},                 // default `<` comparator
            &Particle::id       // project each Particle to its `id`
        );
        int id = it == particle_infos.end() ? 0 : it->id+1; 

        for (int x = -sphere.radius.x; x < sphere.radius.x; ++x) {
            for (int y = -sphere.radius.y; y < sphere.radius.y; ++y) {
                for (int z = -sphere.radius.z; z < sphere.radius.z; ++z) {

                    vec3 pos = {x * width, y * width, z * width};
                    vec3 pos_sq = pos.mul(pos_sq);

                    // if not in ellipsoid skip
                    if ( pos_sq.x / r.x + pos_sq.y / r.y + pos_sq.z / r.z > 1) continue;
                    
                    ids.push_back(id);
                    Particle p =  {
                        .id = id++,
                        .type = sphere.type,
                        .position = sphere.center + pos,
                        .velocity = sphere.mean_velocity,
                        .mass = sphere.mass,
                        .state = sphere.state,
                    };
                    p.velocity += sphere.thermal_velocity(p);

                    add_particle(p);
                }
            }
        }

        return ids;
    }


    void Environment::validate_inputs() {
        //check for duplicates: first sort, then check adjacent pairs
        std::sort(interactions.begin(), interactions.end(), [](const auto& a, const auto& b) {
            if (a.key_pair != b.key_pair)
                return a.key_pair < b.key_pair;
            return a.pair_contains_types < b.pair_contains_types;
        });

        for (size_t i = 1; i < interactions.size(); ++i) {
            if (interactions[i].key_pair == interactions[i-1].key_pair &&
                interactions[i].pair_contains_types == interactions[i-1].pair_contains_types) {
                    throw std::invalid_argument("Found duplicate forces");
            }
        }

        // check if interaction arguments (ids/types) are valid
        for (auto & interaction : interactions) {
            // make sure interaction pair arguments are in canonical order
            auto [a,b] = interaction.key_pair;
            if (a > b)

            if (interaction.pair_contains_types) {
                auto [type1, type2] = interaction.key_pair;
                if (not usr_particle_types.contains(type1) or not usr_particle_types.contains(type2)) {
                    throw std::invalid_argument("Specified interacting particle types do not exist");
                }
            } else {
                auto [id1, id2] = interaction.key_pair;
                if (id1 == PARTICLE_ID_UNDEFINED or id2 == PARTICLE_ID_UNDEFINED) {
                    throw std::invalid_argument("Cannot have interaction between particles with undefined ids");
                }
                if (not usr_particle_ids.contains(id1) or not usr_particle_ids.contains(id2)) {
                    throw std::invalid_argument("Specified interacting particle IDs do not exist");
                }
                if (id1 == id2) {
                    throw std::invalid_argument("Cannot have interaction of a particle with itself");
                }
            }
        }
    }


    void Environment::map_ids_and_types_to_internal() {
    
        // give partciles with undefined id a valid id 
        impl::ParticleID id = 0;
        for (auto & p: particle_infos) {
            if (p.id != PARTICLE_ID_UNDEFINED) continue;

            while (usr_particle_ids.contains(id)) id++;

            usr_particle_ids.insert(id);
            p.id = id;
        }
        
        // generate mapping for user types to implementation types
        std::vector<ParticleType> type_vector;
        type_vector.reserve(usr_particle_types.size());
        type_vector.insert(type_vector.end(), usr_particle_types.begin(), usr_particle_types.end());

        for (size_t  i = 0; i < type_vector.size(); i++) {
            usr_types_to_impl_types[type_vector[i]] = i;
        }

        // generate mapping for user ids to implementation ids
        std::vector<ParticleType> id_vector;
        id_vector.reserve(usr_particle_ids.size());
        id_vector.insert(id_vector.end(), usr_particle_ids.begin(), usr_particle_ids.end());

        // collect all ids involved in an id-to-id interaction 
        std::unordered_set<ParticleID> interacting_ids;
        for (auto & x : interactions) {
            if (!x.pair_contains_types) {
                interacting_ids.insert(x.key_pair.first);
                interacting_ids.insert(x.key_pair.second);
            }
        }
        
        //swap ids, such that all locally interacting particles have the lowest ids
        // for (size_t front = 0, back = id_vector.size()-1; front < back; back--) {
        //     if (interacting_ids.contains(id_vector[back])) {
        //         // find lowest non interacting id
        //         while (interacting_ids.contains(id_vector[front]) && front < back ) {
        //             front++;
        //         }
        //         if (front >= back) break;

        //        std::swap(id_vector[front], id_vector[back]);
        //        front++; 
        //     } 
        // }

        //swap ids, such that all locally interacting particles have the lowest ids
        std::partition(id_vector.begin(), id_vector.end(), 
            [&](ParticleID id) { return interacting_ids.contains(id); }
        );
        
        // create id map
        for (size_t  i = 0; i < id_vector.size(); i++) {
            usr_ids_to_impl_ids[id_vector[i]] = i;
        }
        
    
        // apply mapping to interactions
        // for (auto & interaction : interactions) {
        //     if (interaction.pair_contains_types) {
        //         auto [type1, type2] = interaction.key_pair;
        //         ParticleType impl_type1 = usr_types_to_impl_types[type1];
        //         ParticleType impl_type2 = usr_types_to_impl_types[type2];
        //         interaction.key_pair = {impl_type1, impl_type2}; 
        //     } else {
        //         auto [id1, id2] = interaction.key_pair;
        //         ParticleType impl_id1 = usr_ids_to_impl_ids[id1];
        //         ParticleType impl_id2 = usr_ids_to_impl_ids[id2];
        //         interaction.key_pair = {impl_id1, impl_id2}; 
        //     }
        // }
    }


    void Environment::build_particles() {
        particle_storage.reserve(particle_infos.size());

        size_t idx = 0;
        for (const auto & p : particle_infos) {
            particle_storage.emplace_back(
                idx++, usr_ids_to_impl_ids[p.id], p.position, p.velocity, p.mass,usr_types_to_impl_types[p.type], p.state
            );
        }
    }

    void Environment::build()
    {
        if (is_built) return;

        validate_inputs();
        map_ids_and_types_to_internal();
        build_particles();

        interaction_manager.build(interactions, usr_types_to_impl_types, usr_ids_to_impl_ids);
        
        // free up memory
        interactions.clear();
        particle_infos.clear();
        usr_particle_ids.clear();
        usr_particle_types.clear();
        usr_types_to_impl_types.clear();
        usr_particle_ids.clear();
        usr_ids_to_impl_ids.clear();

        is_built = true;
    }

    impl::ParticleIterator Environment::particles(ParticleState state) {
         return impl::ParticleIterator(particle_storage, state);
    }

    const std::vector<impl::Particle> &Environment::export_particles() {
        return particle_storage;
    }
} // namespace april::env

