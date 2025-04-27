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
		// partition interactions into type and id based interactions
		const auto it = std::partition(interaction_infos.begin(), interaction_infos.end(),
		                               [](const InteractionInfo& info) {
			                               return info.pair_contains_types;
		                               });
		
		std::vector<InteractionInfo> type_interaction_infos {
			std::make_move_iterator(interaction_infos.begin()),
			std::make_move_iterator(it)
		};

		std::vector<InteractionInfo> id_interaction_infos {
			std::make_move_iterator(it),
			std::make_move_iterator(interaction_infos.end())
		};

		// apply types mapping 
		std::vector<std::pair<ParticleType, ParticleType>> type_interaction_keys;
		std::vector<ForcePtr> type_interaction_forces;
		
		type_interaction_keys.reserve(type_interaction_infos.size());
		type_interaction_forces.reserve(type_interaction_infos.size());

		std::unordered_set<ParticleType> particle_types; // needed for mixing forces
		for (auto & x : type_interaction_infos) {
			ParticleType a = usr_types_to_impl_types.at(x.key_pair.first);
			ParticleType b = usr_types_to_impl_types.at(x.key_pair.second);
			
			type_interaction_keys.emplace_back(a,b);
			type_interaction_forces.push_back(std::move(x.force));
			
			particle_types.insert(a);
			particle_types.insert(b);
		}

		// mix type forces
		const size_t n_types = particle_types.size();
		const size_t map_size = 0.5 * (n_types*n_types+n_types); // small gauss sum formula for size of triangle matrix
		std::vector<std::pair<ParticleType, ParticleType>> type_keys_map (map_size); 
		std::vector<ForcePtr> type_forces_map;
		type_forces_map.resize(map_size);

		auto index = [=](const size_t a, const size_t b) {
			return (n_types * a) - (a*(a+1))/2 + (b - a);
		};

		for (size_t i = 0; i < type_interaction_keys.size(); i++) {
			auto [a,b] = type_interaction_keys[i];  //  a <= b implicitly enforced
			const size_t idx = index(a,b); // indexing an upper right triangle
			type_keys_map[idx] = type_interaction_keys[i];
			type_forces_map[idx] = std::move(type_interaction_forces[i]);
		} 
		// ? what if there are particles with types not in any interaction? currently evaulate on force manager will crash
		for (size_t a = 0; a < n_types; a++) {
			for (size_t b = a; b < n_types; b++) {
				const size_t idx =  index(a,b);
				if (type_forces_map[idx] == nullptr) {
					const Force * f1 = type_forces_map[index(a,a)].get();
					const Force * f2 = type_forces_map[index(b,b)].get();

					std::unique_ptr<Force> mixed = f1->mix(f2);
					type_forces_map[index(a,b)] = f1->mix(f2);
					type_keys_map[index(a,b)] = {a,b};
				}
			}
		}

		// apply id mapping
		std::vector<std::pair<ParticleID, ParticleID>> id_interaction_keys;
		std::vector<ForcePtr> id_interaction_forces;
		
		id_interaction_keys.reserve(id_interaction_infos.size());
		id_interaction_forces.reserve(id_interaction_infos.size());

		for (auto & x : id_interaction_infos) {
			ParticleID a = usr_ids_to_impl_ids.at(x.key_pair.first);
			ParticleID b = usr_ids_to_impl_ids.at(x.key_pair.second);
			
			id_interaction_keys.emplace_back(a,b);
			id_interaction_forces.push_back(std::move(x.force));
		}
		
		// build force maps
		inter_type_forces.build(type_keys_map, std::move(type_forces_map));
		intra_particle_forces.build(id_interaction_keys, std::move(id_interaction_forces));
	}


	vec3 InteractionManager::evaluate(const Particle& p1, const Particle& p2, const vec3& distance) const {
		const Force * force_fn = inter_type_forces.get(p1.type, p2.type);
		vec3 force = (*force_fn)(p1, p2, distance); // we always expect a valid function pointer

		if (p1.id < intra_particle_forces.key_size() && p2.id < intra_particle_forces.key_size()) {
			if (const Force * id_force_fn = inter_type_forces.get(p1.type, p2.type)) {
				force += (*id_force_fn)(p1, p2, distance);
			}
		}

		return force;
	}
} 