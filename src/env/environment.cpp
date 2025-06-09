#include "april/env/environment.h"
#include <algorithm>
#include <stdexcept>
#include <ranges>
#include <unordered_set>
#include <format>

#include "april/containers/direct_sum.h"


namespace april::env {
    const auto VEC3_UNDEF = vec3(std::numeric_limits<double>::infinity());

    const vec3 Environment::EXTENT_AUTO = vec3(std::numeric_limits<double>::min());
    const vec3 Environment::ORIGIN_AUTO = vec3(std::numeric_limits<double>::min());

    Environment::Environment():
    is_built(false),
    particle_bbox_extent(0),
    particle_bbox_origin(VEC3_UNDEF),
    container(std::make_unique<core::DirectSum>()),
    interaction_manager()
    {
    }

    void Environment::add_particle(const vec3& position, const vec3& velocity, const double mass, const ParticleType type, const ParticleID id) {
        add_particle(Particle{
            .id = id,
            .type = type,
            .position =  position,
            .velocity = velocity,
            .mass = mass,
            .state = ParticleState::ALIVE
        });
    }

    void Environment::add_particle(const Particle &particle) {
        if (is_built) {
            throw std::logic_error("cannot add particle. environment has already been built.");
        }

        if (particle.id != PARTICLE_ID_UNDEFINED && usr_particle_ids.contains(particle.id)) {
            throw std::invalid_argument("specified id is no unique");
        }

        particle_infos.push_back(particle);   

        // collect id and type information
        if (particle.id != PARTICLE_ID_UNDEFINED) {
            usr_particle_ids.insert(particle.id);
        }
        usr_particle_types.insert(particle.type);

        // adjust bounding box and origin
        if (particle_bbox_origin == VEC3_UNDEF) { //branch is only taken on first particle
            particle_bbox_origin = particle.position;
            return;
        }

        const vec3 &p = particle.position;
        const vec3 old_min = particle_bbox_origin;
        const vec3 old_max = particle_bbox_origin + particle_bbox_extent;

        const vec3 new_min {
            std::min(old_min.x, p.x),
            std::min(old_min.y, p.y),
            std::min(old_min.z, p.z)
        };

        const vec3 new_max {
            std::max(old_max.x, p.x),
            std::max(old_max.y, p.y),
            std::max(old_max.z, p.z)
        };

        particle_bbox_origin = new_min;
        particle_bbox_extent = new_max - new_min;
    }


    void Environment::add_particles(const std::vector<Particle> & particles) {
        for (auto & x: particles) {
            add_particle(x);
        }
    }


    std::vector<ParticleID> Environment::add_particle_cuboid(const ParticleCuboid& cuboid) {
        const uint32_t particle_count = cuboid.particle_count[0] * cuboid.particle_count[1] * cuboid.particle_count[2];
        const double width = cuboid.distance;

        std::vector<ParticleID> ids;
        ids.reserve(particle_count);

        const auto it = std::ranges::max_element(
            particle_infos,
            {},               // default `<` comparator
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
        const double width = sphere.distance;
        const vec3 & r = sphere.radius; 

        std::vector<ParticleID> ids;

        // get the maximum current id
        const auto it = std::ranges::max_element(
            particle_infos,
            {},                 // default `<` comparator
            &Particle::id       // project each Particle to its `id`
        );
        int id = it == particle_infos.end() ? 0 : it->id+1; 

        for (int x = -static_cast<int>(sphere.radius.x/width); x < sphere.radius.x; ++x) {
            for (int y = -static_cast<int>(sphere.radius.y/width); y < sphere.radius.y; ++y) {
                for (int z = -static_cast<int>(sphere.radius.z/width); z < sphere.radius.z; ++z) {

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
        std::ranges::sort(interactions, [](const auto& a, const auto& b) {
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

        // check that each particle type has an interaction
        std::unordered_set<ParticleType> types_without_interaction;
        types_without_interaction.insert(usr_particle_types.begin(), usr_particle_types.end());

        for (const auto & x : interactions) {
            auto [a,b] = x.key_pair;
            if (x.pair_contains_types and a == b and types_without_interaction.contains(a)) {
                types_without_interaction.erase(a);
            }
        }

        if (not types_without_interaction.empty()) {
            std::string msg = "Cannot have particle types without interaction. Types without interaction:";
            for (const auto type : types_without_interaction)
                msg += " " + std::to_string(type);
            throw std::invalid_argument(msg);
        }


        // check if interaction arguments (ids/types) are valid
        for (auto & interaction : interactions) {
            // make sure interaction pair arguments are in canonical order
            if (auto [a,b] = interaction.key_pair; a > b)
                throw std::invalid_argument(
                    "Keys not in canonical order. Pair: (" +
                    std::to_string(a) + ", " + std::to_string(b) + ")"
                );

            if (interaction.pair_contains_types) {
                auto [type1, type2] = interaction.key_pair;
                if (!usr_particle_types.contains(type1) || !usr_particle_types.contains(type2)) {
                    throw std::invalid_argument(
                        "Specified interacting particle types do not exist: (" +
                        std::to_string(type1) + ", " + std::to_string(type2) + ")"
                    );
                }
                if (!usr_particle_types.contains(type1) || !usr_particle_types.contains(type2)) {
                    throw std::invalid_argument(
                        "Specified interacting particle IDs do not exist: (" +
                        std::to_string(type1) + ", " + std::to_string(type2) + ")"
                    );
                }
            } else {
                auto [id1, id2] = interaction.key_pair;
                if (id1 == PARTICLE_ID_UNDEFINED || id2 == PARTICLE_ID_UNDEFINED) {
                    throw std::invalid_argument(
                        "Cannot have interaction between particles with undefined IDs: (" +
                        std::to_string(id1) + ", " + std::to_string(id2) + ")"
                    );
                }
                if (!usr_particle_ids.contains(id1) || !usr_particle_ids.contains(id2)) {
                    throw std::invalid_argument(
                        "Specified interacting particle IDs do not exist: (" +
                        std::to_string(id1) + ", " + std::to_string(id2) + ")"
                    );
                }
                if (id1 == id2) {
                    throw std::invalid_argument(
                        "Cannot have self-interaction of particle ID: " + std::to_string(id1)
                    );
                }
            }
        }

        // check that all particles are contained in the box specified by extent and origin

        // check that extent is larger or equal than particle_bbox_extent
        if (extent != EXTENT_AUTO and (
                std::abs(extent.x) < particle_bbox_extent.x ||
                std::abs(extent.y) < particle_bbox_extent.y ||
                std::abs(extent.z) < particle_bbox_extent.z
            )) {
            throw std::invalid_argument(
                "Specified Environment extent is too small to contain all particles: \n"
                "\tSet Extent (abs): " + vec3{std::abs(extent.x), std::abs(extent.y), std::abs(extent.z)}.to_string() + "\n"
                "\tParticle bounding box size: " + particle_bbox_extent.to_string()
            );
        }

        const vec3 bbox_min = particle_bbox_origin;
        const vec3 bbox_max = particle_bbox_origin + particle_bbox_extent;

        // check that origin is not inside the box specified by (particle_bbox_origin, particle_bbox_extent)
        if (origin != ORIGIN_AUTO and (
            origin.x < bbox_min.x || origin.y < bbox_min.y || origin.z < bbox_min.z ||
            origin.x > bbox_max.x || origin.y > bbox_max.y || origin.z > bbox_max.z
            )) {
            throw std::invalid_argument(
                "Environment origin is not compatible with particle bounding box: \n"
                "\tSet Origin: " + origin.to_string() + "\n"
                "\tParticle bounding box: [" + bbox_min.to_string() + " — " + bbox_max.to_string() + "]"
            );
        }

        //check that (particle_bbox_origin, particle_bbox_extent) is contained in (origin extent)
        if (origin != ORIGIN_AUTO and extent != EXTENT_AUTO) {
            const vec3 env_min {
                std::min(origin.x, origin.x + extent.x),
                std::min(origin.y, origin.y + extent.y),
                std::min(origin.z, origin.z + extent.z)
            };
            const vec3 env_max {
                std::max(origin.x, origin.x + extent.x),
                std::max(origin.y, origin.y + extent.y),
                std::max(origin.z, origin.z + extent.z)
            };

            if (bbox_min.x < env_min.x || bbox_min.y < env_min.y || bbox_min.z < env_min.z ||
                bbox_max.x > env_max.x || bbox_max.y > env_max.y || bbox_max.z > env_max.z)
            {
                throw std::invalid_argument(
                    "Particle cloud lies outside the user‐specified Environment box:\n"
                    "  Env box: [" + env_min.to_string() + " — " + env_max.to_string() + "]\n"
                    "  Particle box: [" + bbox_min .to_string() + " — " + bbox_max .to_string() + "]"
                );
            }
        }
    }


    void Environment::map_ids_and_types_to_internal() {
    
        // give particles with undefined id a valid id
        ParticleID id = 0;
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
            usr_types_to_impl_types[type_vector[i]] = static_cast<impl::ParticleType>(i);
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

        //swap ids, such that all ID-interacting particles have the lowest ids
        std::ranges::partition(id_vector,
                [&](const ParticleID id_) { return interacting_ids.contains(id_); }
        );
        
        // create id map
        for (size_t  i = 0; i < id_vector.size(); i++) {
            usr_ids_to_impl_ids[id_vector[i]] = static_cast<impl::ParticleID>(i);
        }
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

    void Environment::finalize_environment_size() {
        const vec3 bbox_min    = particle_bbox_origin;
        const vec3 bbox_max    = particle_bbox_origin + particle_bbox_extent;
        const vec3 bbox_center = (bbox_min + bbox_max) * 0.5;

        if (extent == EXTENT_AUTO && origin != ORIGIN_AUTO) {
            // User gave origin but no extent:
            // make the box symmetric around `origin` so that bbox_center stays in the middle
            const vec3 opposite_corner = origin + 2 * (bbox_center - origin);
            extent = vec3{
                std::abs(opposite_corner.x - origin.x),
                std::abs(opposite_corner.y - origin.y),
                std::abs(opposite_corner.z - origin.z)
            };
            origin = vec3{
                std::min(origin.x, opposite_corner.x),
                std::min(origin.y, opposite_corner.y),
                std::min(origin.z, opposite_corner.z)
            };
        }

        else if (origin == ORIGIN_AUTO && extent != EXTENT_AUTO) {
            // User gave extent but no origin:
            // center the user‐box on the particle bbox center
            extent = vec3{
                std::abs(extent.x),
                std::abs(extent.y),
                std::abs(extent.z)
            };
            origin = bbox_center - 0.5 * extent;
        }

        else if (origin == ORIGIN_AUTO && extent == EXTENT_AUTO) {
            // Neither origin nor extent given -> default to particle bbox + padding
            extent = particle_bbox_extent * 2.0;        // twice as large in x,y,z
            origin = bbox_center - 0.5 * extent;        // so the box stays centered
        }
    }


    void Environment::set_extent(const vec3& size) {
        if (is_built) {
            throw std::logic_error("cannot set extent. environment has already been built.");
        }

        this->extent = size;
    }

    void Environment::set_origin(const vec3& origin_vec) {
        if (is_built) {
            throw std::logic_error("cannot set origin. environment has already been built.");
        }

        origin = origin_vec;
    }

    void Environment::set_container(std::unique_ptr<core::Container> container_ptr) {
        container = std::move(container_ptr);
    }

    void Environment::build()
    {
        if (is_built) return;

        validate_inputs();
        map_ids_and_types_to_internal();
        build_particles();
        finalize_environment_size();

        interaction_manager.build(interactions, usr_types_to_impl_types, usr_ids_to_impl_ids);
        container->init(&interaction_manager, &particle_storage, extent, origin);
        container->build();

        // free up memory
        interactions.clear();
        particle_infos.clear();
        usr_particle_ids.clear();
        usr_particle_types.clear();
        usr_ids_to_impl_ids.clear();
        usr_types_to_impl_types.clear();

        is_built = true;
    }

    void Environment::update_forces() const {
        container->calculate_forces();
    }

    impl::ParticleIterator Environment::particles(const ParticleState state) {
         return {particle_storage, state};
    }

    const std::vector<impl::Particle> &Environment::export_particles() {
        return particle_storage;
    }
} // namespace april::env

