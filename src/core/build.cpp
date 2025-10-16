#include "april/core/build.h"


namespace april::core::internal {

	using namespace env;
	using namespace container;


	// ---- Domain Validation & Setting ----

	Domain clean_domain(const Domain& domain) {
		// if extent is set to auto then theres nothing to clean
		if (domain.extent == EXTENT_NOT_SET) {
			return domain;
		}
		// if only origin is set to auto then we flip all extent components to positive
		if (domain.origin == ORIGIN_NOT_SET && domain.extent != EXTENT_NOT_SET) {
			const vec3 extent = {
				std::abs(domain.extent.x),
				std::abs(domain.extent.y),
				std::abs(domain.extent.z)
			};

			return {domain.origin, extent };
		}

		// neither ORIGIN_AUTO nor EXTENT_AUTO are set
		const vec3 origin_corner  = domain.origin;
		const vec3 opposite_corner = domain.origin + domain.extent;

		const vec3 min_corner {
			std::min(origin_corner.x, opposite_corner.x),
			std::min(origin_corner.y, opposite_corner.y),
			std::min(origin_corner.x, opposite_corner.z)
		};

		const vec3 max_corner {
			std::max(origin_corner.x, opposite_corner.x),
			std::max(origin_corner.y, opposite_corner.y),
			std::max(origin_corner.x, opposite_corner.z)
		};

		return { min_corner, max_corner - min_corner };
	}


	Box calculate_bounding_box(const std::vector<Particle>& particles) {
		Box bbox;
		if (particles.empty()) return bbox;

		bbox.min = bbox.max = particles[0].position;

		for (const auto& p : particles) {
			bbox.min.x = std::min(bbox.min.x, p.position.x);
			bbox.min.y = std::min(bbox.min.y, p.position.y);
			bbox.min.z = std::min(bbox.min.z, p.position.z);

			bbox.max.x = std::max(bbox.max.x, p.position.x);
			bbox.max.y = std::max(bbox.max.y, p.position.y);
			bbox.max.z = std::max(bbox.max.z, p.position.z);
		}

		bbox.extent = bbox.max - bbox.min;
		return bbox;
	}


	void validate_domain_params(const Domain & domain, const Box & bbox) {
		const vec3& extent = domain.extent;
		const vec3& origin = domain.origin;
		// check that all particles are contained in the box specified by extent and origin

		// 1. sanity: extent large enough
        // check that extent is larger or equal than particle_bbox_extent
        if (extent != EXTENT_NOT_SET and (
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

		// 2. sanity: origin not inside particle box
        // check that origin is not inside the box specified by (particle_bbox_origin, particle_bbox_extent)
        if (origin != ORIGIN_NOT_SET and (
            (origin.x > bbox.min.x && origin.x < bbox.max.x) and
            (origin.y > bbox.min.y && origin.y < bbox.max.y) and
            (origin.z > bbox.min.z && origin.z < bbox.max.z))
           ) {
            throw std::invalid_argument(
                "Environment origin is not compatible with particle bounding box: \n"
                "\tSet Origin: " + origin.to_string() + "\n"
                "\tParticle bounding box: [" + bbox.min.to_string() + " — " + bbox.max.to_string() + "]"
            );
        }

		// 3. containment: ensure all particles fit in environment box
        // check that (particle_bbox_origin, particle_bbox_extent) is contained in (origin extent)
        if (origin != ORIGIN_NOT_SET and extent != EXTENT_NOT_SET) {
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

        	if (bbox.min.x < env_min.x || bbox.min.y < env_min.y || bbox.min.z < env_min.z ||
				bbox.max.x > env_max.x || bbox.max.y > env_max.y || bbox.max.z > env_max.z)
        	{
        		throw std::invalid_argument(
					"Particle cloud lies outside the user‐specified Environment box:\n"
					"  Env box: [" + env_min.to_string() + " — " + env_max.to_string() + "]\n"
					"  Particle box: [" + bbox.min.to_string() + " — " + bbox.max.to_string() + "]"
				);
        	}
        }
	}


	Box finalize_environment_domain( const Box & bbox, const Domain & usr_domain, const vec3 & margin_abs, const vec3 & margin_fac)
	{
		const bool origin_auto = (usr_domain.origin == ORIGIN_NOT_SET);
		const bool extent_auto = (usr_domain.extent == EXTENT_NOT_SET);

		if (margin_abs.x < 0 || margin_abs.y < 0 || margin_abs.z < 0) {
			throw std::logic_error("Absolute margin was set to negative on atleast one axis. Got: " + margin_abs.to_string());
		}

		if (margin_fac.x < 0 || margin_fac.y < 0 || margin_fac.z < 0) {
			throw std::logic_error("Margin factor was set to negative on atleast one axis. Got: " + margin_fac.to_string());
		}

		const vec3 effective_margin = {
			std::max( bbox.extent.x * margin_fac.x, margin_abs.x),
			std::max( bbox.extent.y * margin_fac.y, margin_abs.y),
			std::max( bbox.extent.z * margin_fac.z, margin_abs.z)
		};

		const Box required_box (bbox.min-effective_margin, bbox.max+effective_margin);

		// Case 1: fully manual: both user origin & extent are specified. overrides any margin set
		if (!origin_auto && !extent_auto) {
			return Box(usr_domain);
		}

		// Case 2: fully automatic: both user origin & extent not set
		if (origin_auto && extent_auto) {
			return required_box;
		}

		// Case 3: user origin set, user extent not set
		if (!origin_auto && extent_auto) {
			// we know from validation that origin < bbox.min on all axis so we use that as min corner
			// max corner is chosen such that it satisfies margin requirements
			return {usr_domain.origin, required_box.max};
		}

		// Case 4: user origin not set, user extent set
		if (origin_auto && !extent_auto) {
			// center the particle bounding box inside the simulation domain
			const vec3 bbox_center = (bbox.min + bbox.max) * 0.5;
			const vec3 origin = bbox_center - usr_domain.extent / 2;
			return {origin, origin + usr_domain.extent};
		}

		std::unreachable();
	}




	// ---- Particle Validation & Setting ----

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
		 const std::vector<InteractionParams>& interactions,
		 std::unordered_set<ParticleID> & usr_particle_ids,
		 std::unordered_set<ParticleType> & usr_particle_types
		 ) {
		UserToInternalMappings mapping;

		// give particles with undefined id a valid id
		ParticleID id = 0;
		for (Particle & p: particles) {
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
			mapping.usr_types_to_impl_types[type_vector[i]] = static_cast<env::internal::ParticleType>(i);
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
			mapping.usr_ids_to_impl_ids[id_vector[i]] = static_cast<env::internal::ParticleID>(i);
		}

		return mapping;
	}

	std::vector<env::internal::Particle> build_particles(const std::vector<Particle> & particle_infos, const UserToInternalMappings& mapping) {
		std::vector<env::internal::Particle> particles;
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
