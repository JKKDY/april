#pragma once

#include "april/env/particle.h"
#include "april/env/environment.h"
#include "april/env/interaction.h"
#include "april/algo/container.h"

namespace april::core {

	class System {
	public:
		void update_forces();
		void particles(env::ParticleState state);
		[[nodiscard]] std::vector<env::Particle*> export_particles() const;

	private:
		friend System compile(env::Environment);
		System();

		std::unique_ptr<Container> container;
		env::impl::InteractionManager interaction_manager;
	};


	struct UserToInternalMappings  {
		std::unordered_map<env::ParticleType, env::impl::ParticleType> usr_types_to_impl_types;
		std::unordered_map<env::ParticleID, env::impl::ParticleID> usr_ids_to_impl_ids;
	};


	template<IsContainerDeclaration C> System compile(
		const env::Environment & environment,
		const C & container,
		UserToInternalMappings * particle_mappings = nullptr
	);

}