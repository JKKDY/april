#pragma once
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

		void init(InteractionManager * manager, std::vector<Particle> & particles, const vec3 & extent_vec, const vec3 & origin_vec) {
			interaction_manager = manager;
			particles = particles;
			extent = extent_vec;
			origin = origin_vec;
		}
		virtual void build() = 0;
		virtual void calculate_forces() = 0;
		virtual void update_particle(const Particle &) {}

	protected:
		InteractionManager * interaction_manager{};
		std::vector<Particle> particles{};

		vec3 extent;
		vec3 origin;
	};

	template<typename C> concept IsContainer = std::derived_from<C, Container>;

	template<typename C> concept IsContainerDeclaration = requires {
		typename C::Container;
	}  && IsContainer<typename C::Container>;


} // namespace april::core