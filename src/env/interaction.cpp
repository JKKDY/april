#include "april/env/interaction.h"
#include "april/utils/debug.h"

namespace april::env::impl
{
	// void InteractionManager::register_interaction(const ForceType& force, ParticleTypePair particle_types)
	// {
	// 	AP_ASSERT(particle_types.first < type_count && particle_types.second < type_count,
	// 		"Particle types out of bounds");
	// 	AP_ASSERT(inter_type_forces.size() == type_count * type_count,
	// 		"InteractionManager not initialized properly");

	// 	size_t idx = particle_types.first * type_count + particle_types.second;
	// 	size_t idx2 = particle_types.second * type_count + particle_types.first;

	// 	if (!inter_type_forces[idx]) {
	// 		inter_type_forces[idx] = create_force(force);
	// 	}
	// 	else {
	// 		inter_type_forces[idx] = mix_forces(std::move(inter_type_forces[idx]), force);
	// 	}

	// 	// Create a copy for the symmetric pair
	// 	inter_type_forces[idx2] = std::make_unique<ForceFunctor>(*inter_type_forces[idx]);
	// }




	void InteractionManager::build(const std::vector<Interaction> & interactions) {
		
	}

	vec3 InteractionManager::evaluate(const Particle& p1, const Particle& p2, const vec3& distance) const {
		return vec3();
	}
	


	std::unique_ptr<Force> InteractionManager::mix_forces(ForcePtr force1, ForcePtr force2)
	{
		return std::unique_ptr<Force>();
	}
} 