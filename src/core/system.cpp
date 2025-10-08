#include "april/core/system.h"


namespace april::core::impl {

	using namespace env;
	using namespace container;

	Domain calculate_bounding_box(const std::vector<Particle> & particles) {
		Domain bbox;
		if (particles.empty()) return bbox;

		bbox.origin = particles[0].position;
		for (const auto & particle  : particles) {
			const vec3 &p = particle.position;
			const vec3 old_min = bbox.origin;
			const vec3 old_max = bbox.origin + bbox.extent;

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

			bbox.origin = new_min;
			bbox.extent = new_max - new_min;
		}

		return  bbox;
	}

	void validate_domain_params(const Domain & domain, const Domain & bbox) {
		const vec3& extent = domain.extent;
		const vec3& origin = domain.origin;
		// check that all particles are contained in the box specified by extent and origin
        // first check that extent is larger or equal than particle_bbox_extent
        if (extent != EXTENT_AUTO and (
                std::abs(extent.x) < bbox.extent.x ||
                std::abs(extent.y) < bbox.extent.y ||
                std::abs(extent.z) < bbox.extent.z
            )) {
            throw std::invalid_argument(
                "Specified Environment extent is too small to contain all particles: \n"
                "\tSet Extent (abs): " + vec3{std::abs(extent.x), std::abs(extent.y), std::abs(extent.z)}.to_string() + "\n"
                "\tParticle bounding box size: " + bbox.extent.to_string()
            );
        }

        const vec3 bbox_min = bbox.origin;
        const vec3 bbox_max = bbox.origin + bbox.extent;

        // check that origin is not inside the box specified by (particle_bbox_origin, particle_bbox_extent)
        if (origin != ORIGIN_AUTO and (
            (origin.x > bbox_min.x && origin.x < bbox_max.x) and
            (origin.y > bbox_min.y && origin.y < bbox_max.y) and
            (origin.z > bbox_min.z && origin.z < bbox_max.z))
           ) {
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

	void validate_particle_params(
		const std::vector<Particle> & particles,
		std::vector<InteractionParams> interactions,
		const std::unordered_set<ParticleID> & usr_particle_ids,
		const std::unordered_set<ParticleType> & usr_particle_types)
	{

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

        	// check if types are valid
            if (interaction.pair_contains_types) {
                auto [type1, type2] = interaction.key_pair;
                if (!usr_particle_types.contains(type1) || !usr_particle_types.contains(type2)) {
                    throw std::invalid_argument(
                        "Specified interacting particle types do not exist: (" +
                        std::to_string(type1) + ", " + std::to_string(type2) + ")"
                    );
                }
            } else { // check if ids are valid
                auto [id1, id2] = interaction.key_pair;
                if (id1 == PARTICLE_ID_DONT_CARE || id2 == PARTICLE_ID_DONT_CARE) {
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
                        "Cannot have self-interaction of particle ID: " +
                        std::to_string(id1)
                    );
                }
            }
        }

		// check for positive particle masses
		for (auto & p : particles) {
			if (p.mass <= 0) {
				throw std::invalid_argument(
				   "Particles cannot have negative masses. Particle with ID " +
				   std::to_string(p.id) + " has mass " + std::to_string(p.mass)
				);
			}
		}
	}

	UserToInternalMappings map_ids_and_types_to_internal(
		 std::vector<Particle> & particles,
		 std::vector<InteractionParams> interactions,
		 std::unordered_set<ParticleID> & usr_particle_ids,
		 std::unordered_set<ParticleType> & usr_particle_types
		 ) {
		UserToInternalMappings mapping;

		// give particles with undefined id a valid id
		ParticleID id = 0;
		for (auto & p: particles) {
			if (p.id != PARTICLE_ID_DONT_CARE) continue;

			while (usr_particle_ids.contains(id)) id++;

			usr_particle_ids.insert(id);
			p.id = id;
		}

		// generate mapping for user types to implementation types
		std::vector<ParticleType> type_vector;
		// copy user types
		type_vector.reserve(usr_particle_types.size());
		type_vector.insert(type_vector.end(), usr_particle_types.begin(), usr_particle_types.end());

		// create types map
		for (size_t  i = 0; i < type_vector.size(); i++) {
			mapping.usr_types_to_impl_types[type_vector[i]] = static_cast<env::impl::ParticleType>(i);
		}

		// generate mapping for user ids to implementation ids
		std::vector<ParticleType> id_vector;
		// copy user ids
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

		//swap ids, such that all ID-interacting particles have the lowest implementation ids
		std::ranges::partition(id_vector,
				[&](const ParticleID id_) { return interacting_ids.contains(id_); }
		);

		// create id map
		for (size_t  i = 0; i < id_vector.size(); i++) {
			mapping.usr_ids_to_impl_ids[id_vector[i]] = static_cast<env::impl::ParticleID>(i);
		}

		return mapping;
	}

	Domain finalize_environment_domain(
		const Domain & bbox,
		const Domain & usr_domain
		)
	{
		Domain domain = {usr_domain.extent, usr_domain.origin};
		const vec3& extent = usr_domain.extent;
		const vec3& origin = usr_domain.origin;

		const vec3 bbox_min    = bbox.origin;
		const vec3 bbox_max    = bbox.origin + bbox.extent;
		const vec3 bbox_center = (bbox_min + bbox_max) * 0.5;

		if (extent == EXTENT_AUTO && origin != ORIGIN_AUTO) {
			// User gave origin but no extent:
			// make the box symmetric around `origin` so that bbox_center stays in the middle
			const vec3 opposite_corner = origin + 2 * (bbox_center - origin);
			domain.extent = vec3{
				std::abs(opposite_corner.x - origin.x),
				std::abs(opposite_corner.y - origin.y),
				std::abs(opposite_corner.z - origin.z)
			};
			domain.origin = vec3{
				std::min(origin.x, opposite_corner.x),
				std::min(origin.y, opposite_corner.y),
				std::min(origin.z, opposite_corner.z)
			};

		} else if (origin == ORIGIN_AUTO && extent != EXTENT_AUTO) {
			// User gave extent but no origin:
			// center the user‐box on the particle bbox center
			domain.extent = vec3{
				std::abs(extent.x),
				std::abs(extent.y),
				std::abs(extent.z)
			};
			domain.origin = bbox_center - 0.5 * extent;

		} else if (origin == ORIGIN_AUTO && extent == EXTENT_AUTO) {
			// Neither origin nor extent given -> default to particle bbox + padding
			domain.extent = bbox.extent * 2.0;					// twice as large as bounding box
			domain.origin = bbox_center - 0.5 * domain.extent;  // domain is centered
		}

		AP_ASSERT(domain.extent.x >= bbox.extent.x &&
			domain.extent.y >= bbox.extent.y &&
			domain.extent.z >= bbox.extent.z, "Domain extent must be larger than bounding box");
		return domain;
	}

	std::vector<env::impl::Particle> build_particles(const std::vector<Particle> & particle_infos, const UserToInternalMappings& mapping) {
		std::vector<env::impl::Particle> particles;
		particles.reserve(particle_infos.size());

		for (const auto & p : particle_infos) {
			particles.emplace_back(
				mapping.usr_ids_to_impl_ids.at(p.id),
				p.position,
				p.velocity,
				p.mass,
				mapping.usr_types_to_impl_types.at(p.type),
				p.state,
				vec3{},
				vec3{},
				p.position
			);
		}

		return particles;
	}
}
