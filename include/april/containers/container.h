#pragma once
#include <type_traits>
#include <vector>

#include "april/common.h"

namespace april::env::impl {
	struct Particle;
	class InteractionManager;
}


namespace april::core {

	class Container {
		using InteractionManager = env::impl::InteractionManager;
		using Particle = env::impl::Particle;
	public:
		explicit Container() = default;
		virtual ~Container() = default;

		void init(InteractionManager * interaction_manager, std::vector<Particle> * particles, const vec3 & size, const vec3 & origin) {
			this->interaction_manager = interaction_manager;
			this->particles = particles;
			this->extent = size;
			this->origin = origin;
		}
		virtual void build() = 0;
		virtual void calculate_forces() = 0;
		virtual void update_particle(const Particle &) {}

	protected:
		InteractionManager * interaction_manager{};
		std::vector<Particle> * particles{};

		vec3 extent;
		vec3 origin;
	};

	template<typename T> concept IsContainer = std::is_base_of_v<Container, T>;

} // namespace april::core