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

		void init(InteractionManager * manager, std::vector<Particle> * particles_ptr, const vec3 & size, const vec3 & origin_vec) {
			interaction_manager = manager;
			particles = particles_ptr;
			extent = size;
			origin = origin_vec;
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