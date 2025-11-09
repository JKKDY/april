#pragma once

#include <cxxabi.h>
#include <typeinfo>
#include <algorithm>
#include <unordered_set>
#include <memory>

#include "april/common.h"
#include "april/forces/force.h"
#include "april/particle/descriptors.h"

namespace april::core::internal {
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

	inline void validate_types(
		const std::unordered_set<env::ParticleType>& user_types,
		const std::vector<std::pair<env::ParticleType, env::ParticleType>>& type_pairs
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
		std::unordered_set<env::ParticleType> types_without_interaction;
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


	inline void validate_ids(
		const std::unordered_set<env::ParticleID>& user_ids,
		const std::vector<std::pair<env::ParticleID, env::ParticleID>>& id_pairs
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

	inline void validate_particles(const std::vector<env::Particle>& particles) {
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
	template<force::internal::IsForceVariant FV>
	auto extract_interaction_parameters(
		const std::vector<force::internal::TypeInteraction<FV>> & type_interactions,
		const std::vector<force::internal::IdInteraction<FV>> & id_interaction)
	{
		std::vector<std::pair<env::ParticleType, env::ParticleType>> type_pairs(type_interactions.size());
		std::vector<std::pair<env::ParticleID, env::ParticleID>> id_pairs(id_interaction.size());

		for (size_t i = 0; i < type_interactions.size(); i++)
			type_pairs[i] = {type_interactions[i].type1, type_interactions[i].type2};

		for (size_t i = 0; i < id_interaction.size(); i++)
			id_pairs[i] = {id_interaction[i].id1, id_interaction[i].id2};

		return std::pair {type_pairs, id_pairs};
	}

	inline void assign_missing_particle_ids(
		std::vector<env::Particle>& particles,
		std::unordered_set<env::ParticleID>& user_ids
	) {
		// give particles with undefined id a valid id
		env::ParticleID id = 0;
		for (env::Particle & p: particles) {
			if (p.id.has_value()) continue;

			// find the next free id
			while (user_ids.contains(id)) id++;

			p.id = id;
			user_ids.insert(id);
		}
	}

	inline std::unordered_map<env::ParticleType, env::ParticleType> create_type_map(
		const std::unordered_set<env::ParticleType>& user_types
	) {
		std::unordered_map<env::ParticleType, env::ParticleType> map;
		std::vector<env::ParticleType> type_vector;
		// copy user types
		type_vector.reserve(user_types.size());
		type_vector.insert(type_vector.end(), user_types.begin(), user_types.end());

		// create types map
		for (size_t  i = 0; i < type_vector.size(); i++) {
			map[type_vector[i]] = static_cast<env::ParticleType>(i);
		}

		return map;
	}


	inline std::unordered_map<env::ParticleID, env::ParticleID> create_id_map(
		const std::unordered_set<env::ParticleID>& user_ids,
		const std::vector<std::pair<env::ParticleID, env::ParticleID>>& id_pairs
	) {
		std::unordered_map<env::ParticleID, env::ParticleID> map;
		std::vector<env::ParticleType> id_vector;

		// copy user ids
		id_vector.reserve(user_ids.size());
		id_vector.insert(id_vector.end(), user_ids.begin(), user_ids.end());

		// collect all ids involved in an id-to-id interaction
		std::unordered_set<env::ParticleID> interacting_ids;
		for (const auto & [id1, id2] : id_pairs) {
			interacting_ids.insert(id1);
			interacting_ids.insert(id2);
		}

		//swap ids, such that all ID-interacting particles have the lowest implementation ids
		std::ranges::partition(id_vector,
			[&](const env::ParticleID id_) { return interacting_ids.contains(id_); }
		);

		// create id map
		for (size_t  i = 0; i < id_vector.size(); i++) {
			map[id_vector[i]] = static_cast<env::ParticleID>(i);
		}

		return map;
	}

	inline std::pair<std::unordered_map<env::ParticleType, env::ParticleType>,
	   std::unordered_map<env::ParticleID, env::ParticleID>>
	create_particle_mappings(
		const std::vector<env::Particle>& particles,
		const std::unordered_set<env::ParticleType>& user_types,
		const std::unordered_set<env::ParticleID>& user_ids,
		const std::vector<std::pair<env::ParticleType, env::ParticleType>>& type_pairs,
		const std::vector<std::pair<env::ParticleID, env::ParticleID>>& id_pairs
	) {
		validate_types(user_types, type_pairs);
		validate_ids(user_ids, id_pairs);
		validate_particles(particles);

		auto type_map = create_type_map(user_types);
		auto id_map = create_id_map(user_ids, id_pairs);

		return std::pair{type_map, id_map};
	}

	template<typename T>
		std::string demangled_type_name() {
		int status = 0;
		const std::unique_ptr<char, void(*)(void*)> demangled{
			abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status),
			std::free
		};
		return (status == 0) ? demangled.get() : typeid(T).name();
	}

	template <env::IsUserData UserData>
	std::vector<env::internal::ParticleRecord<UserData>> build_particles(
		const std::vector<env::Particle>& particle_infos,
		const std::unordered_map<env::ParticleType, env::ParticleType> & type_map,
		const std::unordered_map<env::ParticleID, env::ParticleID> & id_map
	) {
		std::vector<env::internal::ParticleRecord<UserData>> particles;
		particles.reserve(particle_infos.size());

		for (const auto & p : particle_infos) {
			AP_ASSERT(p.id.has_value(), "particle id not set during build phase");

			env::internal::ParticleRecord<UserData> particle;
			particle.id = id_map.at(p.id.value());
			particle.type = type_map.at(p.type);
			particle.mass = p.mass;
			particle.state = p.state;
			particle.position = p.position;
			particle.velocity = p.velocity;
			particle.force = p.force.value_or(vec3{});
			particle.old_force = p.old_force.value_or(vec3{});
			particle.old_position = p.old_position.value_or(vec3{});
			if constexpr (std::is_same_v<UserData, env::NoUserData>) {
				particle.user_data = env::NoUserData();
			} else {
				AP_ASSERT((std::any_cast<UserData>(&p.user_data) != nullptr), "user data particle with id "
					+ std::to_string(p.id.value())
					+ " is not of expected type " + demangled_type_name<UserData>()
					+ " but has (mangled) type " + p.user_data.type().name());
				particle.user_data = std::any_cast<UserData>(p.user_data);
			}

			particles.push_back(particle);
		}

		return particles;
	}
}