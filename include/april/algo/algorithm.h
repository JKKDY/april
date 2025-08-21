#pragma once
#include <vector>

#include "april/common.h"
#include "april/env/environment.h"


namespace april::env::impl {
	struct Particle;
	class InteractionManager;
}


namespace april::algo {

	template<typename Config> class Algorithm {
	protected:
		using InteractionManager = env::impl::InteractionManager;
		using Domain = env::Domain;
		using Particle = env::impl::Particle;
	public:
		using CFG = Config;
		Algorithm(Config config, InteractionManager & interactions, const Domain& domain)
			: interactions(interactions), domain(domain), cfg(std::move(config)) {}

		virtual ~Algorithm() = default;
		virtual void build(const std::vector<Particle>& particles) = 0;
		virtual void calculate_forces() = 0;
		virtual void update_particle(const Particle &) {}
	protected:
		InteractionManager & interactions;
		Domain domain;
		Config cfg;
	};

	template<typename A> concept IsAlgo =
		std::derived_from<A, Algorithm<typename A::CFG>> && requires {
	};

	template<typename A> concept IsAlgoDecl = requires {
		typename A::impl;
	}  && IsAlgo<typename A::impl>;


} // namespace april::core