#include "april/core/build.h"


namespace april::core::internal {

	using namespace env;
	using namespace container;


	// ---- Domain Validation & Setting ----

	// Box normalize_domain(const Domain& domain) {
	// 	// if extent is set to auto then there's nothing for us to normalize
	// 	if (domain.extent == EXTENT_NOT_SET) {
	// 		return Box(domain);
	// 	}
	// 	// if only origin is set to auto then we flip all extent components to positive
	// 	if (domain.origin == ORIGIN_NOT_SET && domain.extent != EXTENT_NOT_SET) {
	// 		const vec3 extent = {
	// 			std::abs(domain.extent.x),
	// 			std::abs(domain.extent.y),
	// 			std::abs(domain.extent.z)
	// 		};
	//
	// 		return Box({domain.origin, extent });
	// 	}
	//
	// 	// neither ORIGIN_AUTO nor EXTENT_AUTO are set
	// 	const vec3 origin_corner  = domain.origin;
	// 	const vec3 opposite_corner = domain.origin + domain.extent;
	//
	// 	const vec3 min_corner {
	// 		std::min(origin_corner.x, opposite_corner.x),
	// 		std::min(origin_corner.y, opposite_corner.y),
	// 		std::min(origin_corner.x, opposite_corner.z)
	// 	};
	//
	// 	const vec3 max_corner {
	// 		std::max(origin_corner.x, opposite_corner.x),
	// 		std::max(origin_corner.y, opposite_corner.y),
	// 		std::max(origin_corner.z, opposite_corner.z)
	// 	};
	//
	// 	return {min_corner, max_corner};
	// }


	Box particle_bounding_box(const std::vector<Particle>& particles) {
		if (particles.empty()) return {};

		vec3 min, max = particles[0].position;

		for (const auto& p : particles) {
			min.x = std::min(min.x, p.position.x);
			min.y = std::min(min.y, p.position.y);
			min.z = std::min(min.z, p.position.z);

			max.x = std::max(max.x, p.position.x);
			max.y = std::max(max.y, p.position.y);
			max.z = std::max(max.z, p.position.z);
		}

		return {min, max};
	}


	Box calculate_simulation_box(
		const Domain& desired_domain,
		const Box& required_box,
		const Box& particle_bbox
	) {
		// Case 1: fully manual: both user origin & extent are specified. overrides any margin set
		if (desired_domain.origin.has_value() && desired_domain.extent.has_value()) {
			return Box::from_domain(desired_domain);
		}
		// Case 2: fully automatic: both user origin & extent not set
		if (!desired_domain.origin.has_value() && !desired_domain.extent.has_value()) {
			return required_box;
		}
		// Case 3: user origin set, user extent not set
		if (desired_domain.origin.has_value() && !desired_domain.extent.has_value()) {
			// we know from validation that origin < bbox.min on all axis so we use that as min corner
			// max corner is chosen such that it satisfies margin requirements
			return {desired_domain.origin.value(), required_box.max};
		}
		// Case 4: user origin not set, user extent set
		if (!desired_domain.origin.has_value() && desired_domain.extent.has_value()) {
			// center the particle bounding box inside the simulation domain
			const vec3 bbox_center = (particle_bbox.min + particle_bbox.max) * 0.5;
			const vec3 origin = bbox_center - desired_domain.extent.value() / 2;
			return {origin, origin + desired_domain.extent.value()};
		}
		std::unreachable();
	}


	void verify_domain_consistency(const Box & domain_box, const Box & particle_bbox) {
		// check that extent is larger or equal than particle_bbox_extent
		if (domain_box.extent.x < particle_bbox.extent.x || // extent is always >= 0
			domain_box.extent.y < particle_bbox.extent.y ||
			domain_box.extent.z < particle_bbox.extent.z
		) {
			throw std::invalid_argument(
				"Specified Environment extent is too small to contain all particles: \n"
				"\tDomain box extent: " + domain_box.extent.to_string() + "\n"
				"\tParticle bounding box extent: " + particle_bbox.extent.to_string()
			);
		}

		// ensure all particles fit in environment box
		// check that (particle_bbox_origin, particle_bbox_extent) is contained in (origin extent)
		if (particle_bbox.contains(domain_box.min) || particle_bbox.contains(domain_box.max)) {
			throw std::invalid_argument(
				"Environment origin is not compatible with particle bounding box: \n"
				"\tDomain Box: [" + domain_box.min.to_string() + " — " + domain_box.max.to_string() + "]" "\n"
				"\tParticle bounding box: [" + particle_bbox.min.to_string() + " — " + particle_bbox.max.to_string() + "]"
			);
		}
	}


	Box determine_simulation_box(
		const Domain& desired_domain,
		const Box& particle_bbox,
		const vec3 & margin_abs,
		const vec3 & margin_fac
	) {
		if (margin_abs.x < 0 || margin_abs.y < 0 || margin_abs.z < 0) {
			throw std::logic_error("Absolute margin was set to negative on at least one axis. Got: " + margin_abs.to_string());
		}

		if (margin_fac.x < 0 || margin_fac.y < 0 || margin_fac.z < 0) {
			throw std::logic_error("Margin factor was set to negative on at least one axis. Got: " + margin_fac.to_string());
		}

		const vec3 effective_margin = {
			std::max( particle_bbox.extent.x * margin_fac.x, margin_abs.x),
			std::max( particle_bbox.extent.y * margin_fac.y, margin_abs.y),
			std::max( particle_bbox.extent.z * margin_fac.z, margin_abs.z)
		};

		const Box required_box (particle_bbox.min-effective_margin, particle_bbox.max+effective_margin);

		const Box simulation_box = calculate_simulation_box(desired_domain, required_box, particle_bbox);
		verify_domain_consistency(simulation_box, particle_bbox);

		return {simulation_box.min, simulation_box.max};
	}






	// ---- Particle Validation & Setting ----

	void validate_particle_params(
		const std::vector<Particle> & particles,
		std::vector<InteractionParameters> interactions,
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

	UserToInternalMappings create_particle_mappings(
		 std::vector<Particle> & particles,
		 const std::vector<InteractionParameters>& interactions,
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
			mapping.user_types_to_impl_types[type_vector[i]] = static_cast<env::internal::ParticleType>(i);
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
			mapping.user_ids_to_impl_ids[id_vector[i]] = static_cast<env::internal::ParticleID>(i);
		}

		return mapping;
	}



	std::vector<env::internal::Particle> build_particles(const std::vector<Particle> & particle_infos, const UserToInternalMappings& mapping) {
		std::vector<env::internal::Particle> particles;
		particles.reserve(particle_infos.size());

		for (const auto & p : particle_infos) {
			particles.emplace_back(
				mapping.user_ids_to_impl_ids.at(p.id),
				p.position,
				p.velocity,
				p.mass,
				mapping.user_types_to_impl_types.at(p.type),
				p.state,
				vec3{},
				vec3{},
				p.position
			);
		}

		return particles;
	}


	container::internal::ContainerFlags set_container_flags(const std::vector<boundary::Topology>& topologies) {
		container::internal::ContainerFlags container_flags = {};
		for (const boundary::Face face : boundary::all_faces) {
			if (topologies[face_to_int(face)].force_wrap) {
				switch (axis_of_face(face)) {
				case 0: container_flags.periodic_x = true; break;
				case 1: container_flags.periodic_y = true; break;
				case 2: container_flags.periodic_z = true; break;
				default: std::unreachable();
				}
			}
		}
		return container_flags;
	}
}
