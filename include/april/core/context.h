#pragma once

#include <vector>
#include <cstddef>
#include "april/env/particle.h"
#include "april/env/domain.h"


namespace april::core {

	/// A type erased facade to access the systems functionality.
	/// This is used by monitors, (force) fields & controllers
	/// this allows objects to access the system API without needing
	/// to be templated on it.
	/// While of course slower, than having classes directly templated
	/// on System<...>, it is worth it for the sake of ergonomics, since
	/// these functions should only be called a few times per iteration step
	class SimulationContext {
	public:
		using ParticleView = env::ParticleView;
		using ParticleRef  = env::ParticleRef;

		using ParticleID   = env::internal::ParticleID;
		using Particle     = env::internal::Particle;

		virtual ~SimulationContext() = default;

		// TODO add size(ParticleState) query

		// ---- Core information ----
		[[nodiscard]] virtual env::Domain domain() const noexcept = 0;
		[[nodiscard]] virtual env::Box box() const noexcept = 0;
		[[nodiscard]] virtual double time() const noexcept = 0;
		[[nodiscard]] virtual size_t step() const noexcept = 0;
		[[nodiscard]] virtual size_t size() const noexcept = 0;
		[[nodiscard]] virtual size_t size(env::ParticleState state) const noexcept = 0;


		// ---- Particle access / modification ----
		[[nodiscard]] virtual std::vector<size_t> collect_indices_in_region(const env::Box& region) const = 0;
		[[nodiscard]] virtual std::vector<size_t> collect_indices_in_region(const env::Domain & region) const = 0;

		virtual void register_particle_movement(ParticleID id) = 0;
		virtual void register_all_particle_movements() = 0;

		[[nodiscard]] virtual ParticleRef get_particle_by_id(ParticleID id) noexcept = 0;
		[[nodiscard]] virtual ParticleView get_particle_by_id(ParticleID id) const noexcept = 0;
		[[nodiscard]] virtual ParticleID id_start() const noexcept = 0;
		[[nodiscard]] virtual ParticleID id_end() const noexcept = 0;

		[[nodiscard]] virtual ParticleRef get_particle_by_index(size_t index) noexcept = 0;
		[[nodiscard]] virtual ParticleView get_particle_by_index(size_t index) const noexcept = 0;
		[[nodiscard]] virtual size_t index_start() const noexcept = 0;
		[[nodiscard]] virtual size_t index_end() const noexcept = 0;
	};

	namespace internal {
		/// Adapter template binding a concrete System<> specialization to SimulationContext.
		template<class System> // not restricted due to circular dependencies
		class SimulationContextImpl final : public SimulationContext {
		public:
			explicit SimulationContextImpl(System& sys) : system(sys) {}

			// ---- Core information ----
			[[nodiscard]] env::Domain domain() const noexcept override {
				return system.domain();
			}

			[[nodiscard]] env::Box box() const noexcept override {
				return system.box();
			}

			[[nodiscard]] double time() const noexcept override {
				return system.time();
			}

			[[nodiscard]] size_t step() const noexcept override {
				return system.step();
			}

			[[nodiscard]] size_t size() const noexcept override {
				return system.size();
			}

			[[nodiscard]] size_t size(env::ParticleState) const noexcept override {
				return system.size();
			}



			// ---- Particle modification ----
			void register_particle_movement(ParticleID id) override {
				system.register_particle_movement(id);
			}

			void register_all_particle_movements() override {
				system.register_all_particle_movements();
			}


			// ---- Region Query ----
			[[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Box& region) const override {
				return system.collect_indices_in_region(region);
			}

			[[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Domain & region) const override {
				return system.collect_indices_in_region(region);
			}


			// ---- Particle access by ID ----
			[[nodiscard]] ParticleRef get_particle_by_id(ParticleID id) noexcept override {
				return ParticleRef(system.get_particle_by_id(id));
			}

			[[nodiscard]] ParticleView get_particle_by_id(ParticleID id) const noexcept override {
				return ParticleView(system.get_particle_by_id(id));
			}

			[[nodiscard]] ParticleID id_start() const noexcept override {
				return system.id_start();
			}

			[[nodiscard]] ParticleID id_end() const noexcept override {
				return system.id_end();
			}


			// ---- Particle access by index ----
			[[nodiscard]] ParticleRef get_particle_by_index(size_t index) noexcept override {
				return ParticleRef(system.get_particle_by_index(index));
			}

			[[nodiscard]] ParticleView get_particle_by_index(size_t index) const noexcept override {
				return ParticleView(system.get_particle_by_index(index));
			}

			[[nodiscard]] size_t index_start() const noexcept override {
				return system.index_start();
			}

			[[nodiscard]] size_t index_end() const noexcept override {
				return system.index_end();
			}

		private:
			System& system;
		};
	}

}