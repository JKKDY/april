#pragma once
#include <vector>

#include "april/common.h"
#include "april/env/environment.h"


namespace april::env::impl {
	struct Particle;
	class InteractionManager;
}


namespace april::cont::impl {

	class  IContainer {
	protected:
		using InteractionManager = env::impl::InteractionManager;
		using Domain = env::Domain;
		using Particle = env::impl::Particle;
		using ParticleID = env::impl::ParticleID;
	public:
		IContainer() = default;
		void init(InteractionManager & interaction_mngr, const Domain & dom) {
			interactions = &interaction_mngr;
			domain = dom;
		}
		virtual ~IContainer() = default;

		void dispatch_build(this auto&& self, const std::vector<Particle>& particles) {
	        static_assert(
	            requires { self.build(particles); },
	            "Algorithm subclass must implement: void build(const vector<Particle>&)"
	        );
	        self.build(particles);
	    }

	    void dispatch_calculate_forces(this auto&& self) {
	        static_assert(
	            requires { self.calculate_forces(); },
	            "Algorithm subclass must implement: void calculate_forces()"
	        );
	        self.calculate_forces();
	    }

	    Particle& dispatch_get_particle_by_id(this auto&& self, ParticleID id) {
	        static_assert(
	            requires { { self.get_particle_by_id(id) } -> std::same_as<Particle&>; },
	            "Algorithm subclass must implement: Particle& get_particle_by_id(ParticleID)"
	        );
	        return self.get_particle_by_id(id);
	    }

	    ParticleID dispatch_id_start(this auto&& self) {
	        static_assert(
	            requires { { self.id_start() } -> std::same_as<ParticleID>; },
	            "Algorithm subclass must implement: ParticleID id_start()"
	        );
	        return self.id_start();
	    }

	    ParticleID dispatch_id_end(this auto&& self) {
	        static_assert(
	            requires { { self.id_end() } -> std::same_as<ParticleID>; },
	            "Algorithm subclass must implement: ParticleID id_end()"
	        );
	        return self.id_end();
	    }

	    Particle& dispatch_get_particle_by_index(this auto&& self, size_t index) noexcept {
	        static_assert(
	            requires { { self.get_particle_by_index(index) } -> std::same_as<Particle&>; },
	            "Algorithm subclass must implement: Particle& get_particle_by_index(size_t)"
	        );
	        return self.get_particle_by_index(index);
	    }

	    size_t dispatch_index_start(this auto&& self) {
	        static_assert(
	            requires { { self.index_start() } -> std::same_as<size_t>; },
	            "Algorithm subclass must implement: size_t index_start()"
	        );
	        return self.index_start();
	    }

	    size_t dispatch_index_end(this auto&& self) {
	        static_assert(
	            requires { { self.index_end() } -> std::same_as<size_t>; },
	            "Algorithm subclass must implement: size_t index_end()"
	        );
	        return self.index_end();
	    }

	    size_t dispatch_particle_count(this auto&& self) {
	        static_assert(
	            requires { { self.particle_count() } -> std::same_as<size_t>; },
	            "Algorithm subclass must implement: size_t particle_count()"
	        );
	        return self.particle_count();
	    }
	protected:
		InteractionManager * interactions{};
		Domain domain;
	};


	template<typename Config> class Container : public IContainer {
	public:
		using CFG = Config;
		explicit Container(Config config): cfg(config) {}

	protected:
		Config cfg;
	};


	template<typename A> concept IsContainer =
		std::derived_from<A, Container<typename A::CFG>> && requires {
		typename A::CFG::impl;
	};


	template<typename A> concept IsContDecl = requires {
		typename A::impl;
	}  && IsContainer<typename A::impl>;
} // namespace april::core
