#pragma once
#include <vector>

#include "april/common.h"
#include "april/env/environment.h"


namespace april::env::impl {
	struct Particle;
	class InteractionManager;
}


namespace april::algo::impl {

	class  IAlgorithm {
	protected:
		using InteractionManager = env::impl::InteractionManager;
		using Domain = env::Domain;
		using Particle = env::impl::Particle;
		using ParticleID = env::impl::ParticleID;
	public:
		IAlgorithm() = default;
		void init(InteractionManager & interactions, const Domain & domain) {
			this->interactions = &interactions;
			this->domain = domain;
		}
		virtual ~IAlgorithm() = default;
		virtual void build(const std::vector<env::impl::Particle>& particles) = 0;
		virtual void calculate_forces() = 0;

		virtual Particle & get_particle_by_id(ParticleID id) = 0;
		virtual ParticleID id_start() = 0;
		virtual ParticleID id_end() = 0;

		virtual Particle & get_particle_by_index(size_t index) noexcept = 0;
		virtual size_t index_start() = 0;
		virtual size_t index_end() = 0;

		virtual size_t particle_count() = 0;

	protected:
		InteractionManager * interactions{};
		Domain domain;
	};

	template<typename Config> class Algorithm : public IAlgorithm {
	public:
		using CFG = Config;
		explicit Algorithm(Config config): cfg(config) {}
	protected:
		Config cfg;
	};

	template<typename A> concept IsAlgo =
		std::derived_from<A, Algorithm<typename A::CFG>> && requires {
		typename A::CFG::impl;
	};

	template<typename A> concept IsAlgoDecl = requires {
		typename A::impl;
	}  && IsAlgo<typename A::impl>;
} // namespace april::core
