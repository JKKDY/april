#pragma once

#include "april/containers/container.h"

namespace april::container {

	template<typename Config, typename FT, env::IsUserData U>
	class ContiguousContainer : public Container<Config, FT, U,
	env::internal::ParticleRecordFetcher<U>,
	env::internal::ConstParticleRecordFetcher<U>> {
	public:
		using Base = Container<Config, FT, U,
		env::internal::ParticleRecordFetcher<U>,
		env::internal::ConstParticleRecordFetcher<U>>;
		using typename Base::MutableFetcher;
		using typename Base::ConstFetcher;
		using typename Base::ParticleRecord;
		using Base::Base;

		void build_storage(const std::vector<ParticleRecord>& particles) {
			AP_ASSERT(!is_built, "storage has already been built");

			this->particles = std::vector(particles);
			indices.resize(particles.size());
			for (size_t i = 0; i < particles.size(); i++) {
				const auto id = static_cast<size_t>(this->particles[i].id);
				indices[id] = i;
			}
			this->is_built = true;
		}

		[[nodiscard]] MutableFetcher get_fetcher_by_id(const env::ParticleID id) noexcept{
			// TODO this method needs testing
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			AP_ASSERT(id >= id_start() && id < id_end(), "invalid id. got " + std::to_string(id));
			const size_t index = id_to_index(id);
			return MutableFetcher(particles[index]);
		}

		[[nodiscard]] env::ParticleID id_start() const {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			return 0;
		}

		[[nodiscard]] env::ParticleID id_end() const {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			return static_cast<env::ParticleID>(particles.size());
		}

		[[nodiscard]] size_t id_to_index(const env::ParticleID id) const {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			return indices[static_cast<size_t>(id)];
		}

		// index for non-stable iteration
		[[nodiscard]] MutableFetcher get_fetcher_by_index(size_t index) noexcept {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called beforehand");
			AP_ASSERT(index < particles.size(), "index must be < #particles");
			return MutableFetcher(particles[index]);
		}

		[[nodiscard]] ConstFetcher get_fetcher_by_index(size_t index) const noexcept{
			AP_ASSERT(is_built, "storage was not built. build_storage must be called beforehand");
			AP_ASSERT(index < particles.size(), "index must be < #particles");
			return ConstFetcher(particles[index]);
		}

		[[nodiscard]] size_t index_start() const {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			return 0;
		}

		[[nodiscard]] size_t index_end() const {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			return particles.size();
		}

		[[nodiscard]] size_t particle_count() const {
			return particles.size();
		}

	protected:
		void swap_particles (uint32_t idx1, uint32_t idx2) {
			auto id1 = particles[idx1].id, id2 = particles[idx2].id;
			std::swap(particles[idx1], particles[idx2]);
			std::swap(indices[id1], indices[id2]);
		}

		std::vector<ParticleRecord> particles = {};
		std::vector<uint32_t> indices; // map id to index
	private:
		bool is_built = false;
	};
}