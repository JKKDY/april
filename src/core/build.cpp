#include "april/core/build.h"
#include <algorithm>

namespace april::core::internal {

	using namespace env;
	using namespace container;


	// ---- Domain Validation & Setting ----

	Box particle_bounding_box(const std::vector<Particle>& particles) {
		if (particles.empty()) return {};

		vec3 min = particles[0].position;
		vec3 max = particles[0].position;

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


	void verify_domain_consistency(const Box & simulation_box, const Box & particle_bbox) {
		if (simulation_box.extent.x < 0 || simulation_box.extent.y < 0 || simulation_box.extent.z < 0)
		{
			throw std::logic_error("Simulation domain has negative extent. Got extent " + simulation_box.extent.to_string());
		}

		if (simulation_box.extent.x == 0 && simulation_box.extent.y == 0 && simulation_box.extent.z == 0)
		{
			throw std::logic_error("Simulation domain size is zero. Got extent " + simulation_box.extent.to_string() +
				"\n If you have no particles or they all have the same position, you must specify a domain manually.");
		}

		// check that min corner is outside of particle bbox on all axis
		if (simulation_box.min.x > particle_bbox.min.x ||
			simulation_box.min.y > particle_bbox.min.y ||
			simulation_box.min.z > particle_bbox.min.z
		) {
			throw std::invalid_argument(
				"Specified Environment domain does not contain all particles: \n"
				"\tDomain box min corner: " + simulation_box.min.to_string() + "\n"
				"\tParticle bounding min corner: " + particle_bbox.min.to_string()
			);
		}

		// check that max corner is outside of particle bbox on all axis
		if (simulation_box.max.x < particle_bbox.max.x ||
			simulation_box.max.y < particle_bbox.max.y ||
			simulation_box.max.z < particle_bbox.max.z
		) {
			throw std::invalid_argument(
				"Specified Environment domain does not contain all particles: \n"
				"\tDomain box max corner: " + simulation_box.max.to_string() + "\n"
				"\tParticle bounding max corner: " + particle_bbox.max.to_string()
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



	// ---- Particle Validation ----

	template <typename T>
	void validate_no_duplicates(const std::vector<std::pair<T, T>>& pairs, std::string_view item_name) {
		//check for duplicates: first sort, then check adjacent pairs
		auto sorted_pairs = pairs;
		std::ranges::sort(sorted_pairs);

		// returns the first pair of adjacent, equal elements
		auto it = std::ranges::adjacent_find(sorted_pairs);

		if (it != sorted_pairs.end()) {
			throw std::invalid_argument(std::format(
				"Found duplicate {}: ({}, {})",
				item_name, it->first, it->second
			));
		}
	}

	void validate_types(
		const std::unordered_set<ParticleType>& user_types,
		const std::vector<std::pair<ParticleType, ParticleType>>& type_pairs
	) {
		validate_no_duplicates(type_pairs, "type interaction");

		// check if types are valid i.e. if there are any particles with that type
		for (const auto & [type1, type2] : type_pairs) {
			if (!user_types.contains(type1))
				throw std::invalid_argument(
					"Specified interacting particle type does not exist: " + std::to_string(type1));

			if (!user_types.contains(type2))
				throw std::invalid_argument(
					"Specified interacting particle type does not exist: " + std::to_string(type2));
		}

		// check that each particle type has a self interaction
		std::unordered_set<ParticleType> types_without_interaction;
		types_without_interaction.insert(user_types.begin(), user_types.end());

		for (const auto & [a,b] : type_pairs) {
			if (a == b and types_without_interaction.contains(a)) {
				types_without_interaction.erase(a);
			}
		}

		if (not types_without_interaction.empty()) {
			std::string msg = "Cannot have particle types without interaction. Types without interaction:";
			for (const auto type : types_without_interaction)
				msg += " " + std::to_string(type);
			throw std::invalid_argument(msg);
		}
	}

	void validate_ids(
		const std::unordered_set<ParticleID>& user_ids,
		const std::vector<std::pair<ParticleID, ParticleID>>& id_pairs
	) {
		validate_no_duplicates(id_pairs, "ID interaction");

		// check if interaction arguments (ids/types) are valid
		for (const auto & [id1, id2]  : id_pairs) {

			// check if types are valid
			if (!user_ids.contains(id1) || !user_ids.contains(id2)) {
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

	void validate_particles(const std::vector<Particle>& particles) {
		// check for positive particle masses
		for (auto & p : particles) {
			if (p.mass <= 0) {
				throw std::invalid_argument(
				   "Particles cannot have negative masses. Particle with ID " +
				   std::to_string(p.id.value()) + " has mass " + std::to_string(p.mass)
				);
			}
		}
	}


	// ---- build particles ----

	void assign_missing_particle_ids(
		std::vector<Particle>& particles,
		std::unordered_set<ParticleID>& user_ids
	) {
		// give particles with undefined id a valid id
		ParticleID id = 0;
		for (Particle & p: particles) {
			if (p.id.has_value()) continue;

			// find the next free id
			while (user_ids.contains(id)) id++;

			p.id = id;
			user_ids.insert(id);
		}
	}

	std::unordered_map<ParticleType, ParticleType> create_type_map(
		const std::unordered_set<ParticleType>& user_types
	) {
		std::unordered_map<ParticleType, ParticleType> map;
		std::vector<ParticleType> type_vector;
		// copy user types
		type_vector.reserve(user_types.size());
		type_vector.insert(type_vector.end(), user_types.begin(), user_types.end());

		// create types map
		for (size_t  i = 0; i < type_vector.size(); i++) {
			map[type_vector[i]] = static_cast<ParticleType>(i);
		}

		return map;
	}

	std::unordered_map<ParticleID, ParticleID> create_id_map(
		const std::unordered_set<ParticleID>& user_ids,
		const std::vector<std::pair<ParticleID, ParticleID>>& id_pairs
	) {
		std::unordered_map<ParticleID, ParticleID> map;
		std::vector<ParticleType> id_vector;

		// copy user ids
		id_vector.reserve(user_ids.size());
		id_vector.insert(id_vector.end(), user_ids.begin(), user_ids.end());

		// collect all ids involved in an id-to-id interaction
		std::unordered_set<ParticleID> interacting_ids;
		for (const auto & [id1, id2] : id_pairs) {
			interacting_ids.insert(id1);
			interacting_ids.insert(id2);
		}

		//swap ids, such that all ID-interacting particles have the lowest implementation ids
		std::ranges::partition(id_vector,
			[&](const ParticleID id_) { return interacting_ids.contains(id_); }
		);

		// create id map
		for (size_t  i = 0; i < id_vector.size(); i++) {
			map[id_vector[i]] = static_cast<ParticleID>(i);
		}

		return map;
	}

	std::pair<std::unordered_map<ParticleType, ParticleType>,
			std::unordered_map<ParticleID, ParticleID>>
	create_particle_mappings(
		const std::vector<Particle>& particles,
		const std::unordered_set<ParticleType>& user_types,
		const std::unordered_set<ParticleID>& user_ids,
		const std::vector<std::pair<ParticleType, ParticleType>>& type_pairs,
		const std::vector<std::pair<ParticleID, ParticleID>>& id_pairs
	) {
		validate_types(user_types, type_pairs);
		validate_ids(user_ids, id_pairs);
		validate_particles(particles);

		auto type_map = create_type_map(user_types);
		auto id_map = create_id_map(user_ids, id_pairs);

		return std::pair{type_map, id_map};
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


	void validate_topologies(const std::vector<boundary::Topology> &) {

	}

}
