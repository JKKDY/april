#pragma once

#include <vector>
#include <cstddef>
#include "april/particle/fields.hpp"
#include "april/env/domain.hpp"


namespace april::core {

	template<class System>
	class SystemContext {
	public:

		// -----------------
		// LIFECYCLE & STATE
		// -----------------
		explicit SystemContext(System & sys): system(sys) {}
		[[nodiscard]] double time() const noexcept { return system.time(); }
		[[nodiscard]] size_t step() const noexcept { return system.step(); }
		[[nodiscard]] env::Box box() const noexcept { return system.box(); }


		// ------------------
		// PARTICLE ACCESSORS
		// ------------------
		// INDEX ACCESSORS (fast)
		template<ParticleField M>
		[[nodiscard]] auto at(size_t index) {
			return system.template at<M>(index);
		}

		template<ParticleField M>
		[[nodiscard]] auto view(size_t index) const {
			return system.template view<M>(index);
		}

		template<ParticleField M>
		[[nodiscard]] auto restricted_at(size_t index) {
			return system.template restricted_at<M>(index);
		}

		// ID ACCESSORS (stable)
		template<ParticleField M>
		[[nodiscard]] auto at_id(ParticleID id) {
			return system.template at_id<M>(id);
		}

		template<ParticleField M>
		[[nodiscard]] auto view_id(ParticleID id) const {
			return system.template view_id<M>(id);
		}

		template<ParticleField M>
		[[nodiscard]] auto restricted_at_id(ParticleID id) {
			return system.template restricted_at_id<M>(id);
		}


		// -----------
		// ID INDEXING
		// -----------
		[[nodiscard]] ParticleID min_id() const noexcept{
			return system.min_id();
		}

		// get the largest particle id
		[[nodiscard]] ParticleID max_id() const noexcept {
			return system.max_id();
		}

		[[nodiscard]] bool contains(ParticleID id) const {
			return system.contains(id);
		}


		// -------
		// QUERIES
		// -------
		[[nodiscard]] size_t size(ParticleState = ParticleState::ALL) const noexcept {
			return system.size();
		}

		[[nodiscard]] std::vector<size_t> query_region(const env::Box & region) const {
			return system.query_region(region);
		}

		[[nodiscard]] std::vector<size_t> query_region(const env::Domain & region) const {
			return query_region(env::Box(region.min_corner().value(), region.max_corner().value()));
		}


		// --------------
		// FUNCTIONAL OPS
		// --------------
		template<ParticleField M, ExecutionPolicy Policy = ExecutionPolicy::Seq, typename Func>
		void for_each_particle(Func && func, ParticleState state = ParticleState::ALL) {
			system.template for_each_particle<M, Policy, Func>(std::forward<Func>(func), state);
		}

		template<ParticleField M, ExecutionPolicy Policy = ExecutionPolicy::Seq, typename Func>
		void for_each_particle_view(Func && func, ParticleState state = ParticleState::ALL) const {
			system.template for_each_particle_view<M, Policy, Func>(std::forward<Func>(func), state);
		}

		template<typename Func>
		void for_each_interaction_batch(Func && func) {
			system.for_each_interaction_batch(std::forward<Func>(func));
		}

		template<ParticleField M, typename Func>
		void for_each_interaction_pair(Func && func) {
			system.template for_each_interaction_pair<M>(func);
		}


		// -----------------
		// STRUCTURE UPDATES
		// -----------------
		void rebuild_structure() {
			system.rebuild_structure();
		}

		void notify_moved(const std::vector<size_t> & indices) {
			system.notify_moved(indices);
		}

		void notify_moved_id(const std::vector<ParticleID> & ids) {
			system.notify_moved_id(ids);
		}

	private:
		System& system;
	};

}



