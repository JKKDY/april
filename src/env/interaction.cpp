#include "april/env/interaction.h"
#include "april/utils/debug.h"

#include <algorithm>

namespace april::env::impl
{
	
	void InteractionManager::build(std::vector<InteractionInfo> & interactions, 
		const std::unordered_map<env::ParticleType, impl::ParticleType> & usr_types_to_impl_types,
		const std::unordered_map<env::ParticleID, impl::ParticleID> & usr_ids_to_impl_ids
	) {		
		// partition interactions into type and id based interactions
		auto it = std::partition(interactions.begin(), interactions.end(),
                         [](const InteractionInfo& info) {
                             return info.pair_contains_types;
                         }
					);
		
		std::vector<InteractionInfo> type_interaction{
			std::make_move_iterator(interactions.begin()),
			std::make_move_iterator(it)
		};

		std::vector<InteractionInfo> id_interaction{
			std::make_move_iterator(it),
			std::make_move_iterator(interactions.end())
		};

		// apply usr_types_to_impl_types mapping and 
		std::vector<std::pair<ParticleType, ParticleType>> type_keys;
		type_keys.reserve(type_interaction.size());
		std::vector<ForcePtr> type_forces;
		type_forces.reserve(type_interaction.size());

		for (auto & x : type_interaction) {
			ParticleType a = usr_types_to_impl_types.at(x.key_pair.first);
			ParticleType b = usr_types_to_impl_types.at(x.key_pair.second);
			type_keys.emplace_back(a,b);
			type_forces.push_back(std::move(x.force));
		}

		// apply usr_ids_to_impl_ids mapping
		std::vector<std::pair<ParticleID, ParticleID>> id_keys;
		id_keys.reserve(id_interaction.size());
		std::vector<ForcePtr> id_forces;
		id_forces.reserve(id_interaction.size());

		for (auto & x : id_interaction) {
			ParticleID a = usr_ids_to_impl_ids.at(x.key_pair.first);
			ParticleID b = usr_ids_to_impl_ids.at(x.key_pair.second);
			id_keys.emplace_back(a,b);
			id_forces.push_back(std::move(x.force));
		}

		// build force maps
		inter_type_forces.build(type_keys, std::move(type_forces));
		intra_particle_forces.build(id_keys, std::move(id_forces));
	}

	vec3 InteractionManager::evaluate(const Particle& p1, const Particle& p2, const vec3& distance) const {
		Force* force_fn = inter_type_forces.get(p1.type, p2.type);
		vec3 force = (*force_fn)(p1, p2, distance);

		if (p1.id < intra_particle_forces.key_size() && p2.id < intra_particle_forces.key_size()) {
			Force* id_force_fn = inter_type_forces.get(p1.type, p2.type);
			if (id_force_fn) {
				force += (*id_force_fn)(p1, p2, distance);
			}
		}

		return force;
	}
	


	ForcePtr InteractionManager::mix_forces(ForcePtr force1, ForcePtr force2)
	{
		return std::unique_ptr<Force>();
	}
} 