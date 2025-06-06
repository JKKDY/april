#include "april/env/interaction.h"

#include <algorithm>
#include <unordered_set>
#include <memory>
#include <utility>


namespace april::env::impl
{
	
	void InteractionManager::build(std::vector<InteractionInfo> & interaction_infos, 
		const std::unordered_map<env::ParticleType, impl::ParticleType> & usr_types_to_impl_types,
		const std::unordered_map<env::ParticleID, impl::ParticleID> & usr_ids_to_impl_ids
	) {
		// get the max cutoff distance
		const auto max_it = std::ranges::max_element(
			interaction_infos,
			{},
			[](const InteractionInfo& i) { return i.force->cutoff_radius; }
		);
		max_cutoff = (max_it != interaction_infos.end()) ? max_it->force->cutoff_radius : 0.0;

		// partition interactions into type and id based interactions
		const auto it = std::partition(interaction_infos.begin(), interaction_infos.end(),
		                               [](const InteractionInfo& info) {
			                               return info.pair_contains_types;
		                               });

		// contains all interaction infos for particle types
		std::vector<InteractionInfo> type_interaction_infos {
			std::make_move_iterator(interaction_infos.begin()),
			std::make_move_iterator(it)
		};

		// contains all interaction infos for particle (id) pairs
		std::vector<InteractionInfo> id_interaction_infos {
			std::make_move_iterator(it),
			std::make_move_iterator(interaction_infos.end())
		};

		// apply types mapping (user types -> internal types)
		std::vector<std::pair<ParticleType, ParticleType>> type_interaction_keys_mapping;
		std::vector<ForcePtr> type_interaction_forces_mapping;

		type_interaction_keys_mapping.reserve(type_interaction_infos.size());
		type_interaction_forces_mapping.reserve(type_interaction_infos.size());

		std::unordered_set<ParticleType> particle_types; // needed to perform forces mixing
		for (auto & x : type_interaction_infos) {
			// apply mappings
			ParticleType a = usr_types_to_impl_types.at(x.key_pair.first);
			ParticleType b = usr_types_to_impl_types.at(x.key_pair.second);

			type_interaction_keys_mapping.emplace_back(a,b);
			type_interaction_forces_mapping.push_back(std::move(x.force));

			particle_types.insert(a);
			particle_types.insert(b);
		}

		// mix type forces
		const size_t n_types = particle_types.size();
		const size_t map_size = (n_types*n_types+n_types) / 2; // small gauss sum formula for number of all pairs
		std::vector<std::pair<ParticleType, ParticleType>> type_interaction_keys (map_size);
		std::vector<ForcePtr> type_interaction_forces;
		type_interaction_forces.resize(map_size);

		auto index = [=](const size_t a, const size_t b) {  // indexing an upper right triangle
			// first part: base offset, second part: sum of elements not part of the triangle
			return ((n_types + 1) * a) - (a*(a+1))/2 + (b - a);
		};

		// add all user specified interactions between particle types
		for (size_t i = 0; i < type_interaction_keys_mapping.size(); i++) {
			auto [a,b] = type_interaction_keys_mapping[i];  //  a <= b implicitly enforced
			type_interaction_keys[index(a,b)] = type_interaction_keys_mapping[i];
			type_interaction_forces[index(a,b)] = std::move(type_interaction_forces_mapping[i]);
		}

		// add missing particle type interactions by mixing forces
		for (size_t a = 0; a < n_types; a++) {
			for (size_t b = a; b < n_types; b++) {
				const size_t idx =  index(a,b);
				if (type_interaction_forces[idx] == nullptr) {
					const Force * f1 = type_interaction_forces[index(a,a)].get();
					const Force * f2 = type_interaction_forces[index(b,b)].get();

					//? will the asymmetry that we are only calling f1->mix become a problem?
					std::unique_ptr<Force> mixed = f1->mix(f2);
					type_interaction_forces[index(a,b)] = f1->mix(f2);
					type_interaction_keys[index(a,b)] = {static_cast<ParticleType>(a),static_cast<ParticleType>(b)};
				}
			}
		}

		// apply id mapping (user specified ids -> internal ids)
		std::vector<std::pair<ParticleID, ParticleID>> id_interaction_keys;
		std::vector<ForcePtr> id_interaction_forces;

		id_interaction_keys.reserve(id_interaction_infos.size());
		id_interaction_forces.reserve(id_interaction_infos.size());

		for (auto & x : id_interaction_infos) {
			// apply mappings
			ParticleID a = usr_ids_to_impl_ids.at(x.key_pair.first);
			ParticleID b = usr_ids_to_impl_ids.at(x.key_pair.second);

			id_interaction_keys.emplace_back(a,b);
			id_interaction_forces.push_back(std::move(x.force));
		}

		// build force maps
		inter_type_forces.build(type_interaction_keys, std::move(type_interaction_forces));
		intra_particle_forces.build(id_interaction_keys, std::move(id_interaction_forces));
	}

	vec3 InteractionManager::evaluate(const Particle& p1, const Particle& p2) const {
		return evaluate(p1, p2, p2.position - p1.position); // dist vector points from p1 to p2
	}


	vec3 InteractionManager::evaluate(const Particle& p1, const Particle& p2, const vec3& distance) const {
		const Force * force_fn = inter_type_forces.get(p1.type, p2.type);
		vec3 force = (*force_fn)(p1, p2, distance); // we always expect a valid function pointer

		// check if both particles even have any individual interactions defined for them
		if (p1.id < intra_particle_forces.key_size() && p2.id < intra_particle_forces.key_size()) {
			if (const Force * id_force_fn = intra_particle_forces.get(p1.id, p2.id)) {
				force += (*id_force_fn)(p1, p2, distance);
			}
		}

		return force;
	}

	double InteractionManager::get_max_cutoff() const {
		return max_cutoff;
	}
} 