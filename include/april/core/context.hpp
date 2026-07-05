#pragma once

#include <vector>
#include "domain.hpp"
#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"
#include "april/particle/particle_types.hpp"

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
		[[nodiscard]] Box box() const noexcept { return system.box(); }


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


		// ID ACCESSORS (stable)
		template<ParticleField M>
		[[nodiscard]] auto at_id(ParticleID id) {
			return system.template at_id<M>(id);
		}

		template<ParticleField M>
		[[nodiscard]] auto view_id(ParticleID id) const {
			return system.template view_id<M>(id);
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

		[[nodiscard]] bool contains_id(ParticleID id) const {
			return system.contains_id(id);
		}


		// -------
		// QUERIES
		// -------
		[[nodiscard]] size_t size(ParticleState state = ParticleState::ALL) const noexcept {
			return system.size(state);
		}

		[[nodiscard]] std::vector<size_t> query_region(const core::Box & region) const {
			return system.query_region(region);
		}

		[[nodiscard]] std::vector<size_t> query_region(const Domain & region) const {
			return query_region(core::Box(region.min_corner().value(), region.max_corner().value()));
		}


		// --------------
		// FUNCTIONAL OPS
		// --------------
		template<
			ParallelPolicy P = ParallelPolicy::Serial,
			VectorPolicy V = VectorPolicy::Auto,
			exec::IsKernel Kernel>
		void for_each_particle(Kernel && func, ParticleState state = ParticleState::ALL) {
			system.template for_each_particle<P, V, Kernel>(std::forward<Kernel>(func), state);
		}

		template<
			ParallelPolicy P = ParallelPolicy::Serial,
			VectorPolicy V = VectorPolicy::Auto,
			exec::IsKernel Kernel>
		void for_each_particle_view(Kernel && func, ParticleState state = ParticleState::ALL) const {
			system.template for_each_particle_view<P, V, Kernel>(std::forward<Kernel>(func), state);
		}

		template<exec::IsKernel Func>
		void for_each_interaction_batch(Func && func) {
			system.for_each_interaction_batch(std::forward<Func>(func));
		}

		template<exec::IsKernel Func>
		void for_each_interaction_pair(Func && func) {
			system.for_each_interaction_pair(func);
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














