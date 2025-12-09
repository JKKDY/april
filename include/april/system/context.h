#pragma once

#include <vector>
#include <cstddef>
#include "april/particle/fields.h"
#include "april/env/domain.h"


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
		template<env::FieldMask M>
		[[nodiscard]] auto at(size_t index) {
			return system.template at<M>(index);
		}

		template<env::FieldMask M>
		[[nodiscard]] auto view(size_t index) const {
			return system.template view<M>(index);
		}

		template<env::FieldMask M>
		[[nodiscard]] auto restricted_at(size_t index) {
			return system.template restricted_at<M>(index);
		}

		// ID ACCESSORS (stable)
		template<env::FieldMask M>
		[[nodiscard]] auto at_id(env::ParticleID id) {
			return system.template at_id<M>(id);
		}

		template<env::FieldMask M>
		[[nodiscard]] auto view_id(env::ParticleID id) const {
			return system.template view_id<M>(id);
		}

		template<env::FieldMask M>
		[[nodiscard]] auto restricted_at_id(env::ParticleID id) {
			return system.template restricted_at_id<M>(id);
		}


		// -----------
		// ID INDEXING
		// -----------
		[[nodiscard]] env::ParticleID min_id() const noexcept{
			return system.min_id();
		}

		// get the largest particle id
		[[nodiscard]] env::ParticleID max_id() const noexcept {
			return system.max_id();
		}

		[[nodiscard]] bool contains(env::ParticleID id) const {
			return system.contains(id);
		}


		// -------
		// QUERIES
		// -------
		[[nodiscard]] size_t size(env::ParticleState = env::ParticleState::ALL) const noexcept {
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
		template<env::FieldMask M, typename Func, bool parallelize=false>
		void for_each_particle(Func && func) {
			system.template for_each_particle<M, Func, parallelize>(std::forward<Func>(func));
		}

		template<typename Func>
		void for_each_interaction_batch(Func && func) {
			system.for_each_interaction_batch(std::forward<Func>(func));
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

		void notify_moved_id(const std::vector<env::ParticleID> & ids) {
			system.notify_moved_id(ids);
		}

	private:
		System& system;
	};

}