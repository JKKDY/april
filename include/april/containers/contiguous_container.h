#pragma once

#include "april/particle/particle.h"
#include "april/containers/container.h"

namespace april::container {

	namespace internal {
		template<env::IsUserData U>
		struct ParticleRecordFetcher {
			using UserDataT = U;
			using Record = env::internal::ParticleRecord<UserDataT>;

			explicit ParticleRecordFetcher(Record& r) : record(r) {}

			[[nodiscard]] vec3& position()       		const noexcept { return record.position; }
			[[nodiscard]] vec3& velocity()       		const noexcept { return record.velocity; }
			[[nodiscard]] vec3& force()          		const noexcept { return record.force; }
			[[nodiscard]] vec3& old_position()   		const noexcept { return record.old_position; }
			[[nodiscard]] vec3& old_force()      		const noexcept { return record.old_force; }
			[[nodiscard]] double& mass()         		const noexcept { return record.mass; }
			[[nodiscard]] env::ParticleState& state()	const noexcept { return record.state; }
			[[nodiscard]] env::ParticleType& type()		const noexcept { return record.type; }
			[[nodiscard]] env::ParticleID& id()			const noexcept { return record.id; }
			[[nodiscard]] UserDataT& user_data()		const noexcept { return record.user_data; }
		private:
			Record& record;
		};

	}


	template<typename Config, env::IsUserData U>
	class ContiguousContainer : public Container<Config, U, internal::ParticleRecordFetcher<U>>{
	public:
		using Base = Container<Config, U, internal::ParticleRecordFetcher<U>>;
		friend Base;

		using Particle = env::internal::ParticleRecord<U>;
		using Base::Base;
		using typename Base::FetcherT;

		void build_storage(const std::vector<Particle>& particles) {
			AP_ASSERT(!is_built, "storage has already been built");

			this->particles = std::vector(particles);
			id_to_index_map.resize(particles.size());
			for (size_t i = 0; i < particles.size(); i++) {
				const auto id = static_cast<size_t>(particles[i].id);
				id_to_index_map[id] = i;
			}
			this->is_built = true;
		}

		// INDEXING
		[[nodiscard]] size_t id_to_index(const env::ParticleID id) const {
			return id_to_index_map[static_cast<size_t>(id)];
		}

		[[nodiscard]] env::ParticleID min_id() const {
			return 0;
		}

		[[nodiscard]] env::ParticleID max_id() const {
			return static_cast<env::ParticleID>(particles.size());
		}


		// QUERIES
		[[nodiscard]] bool contains(const env::ParticleID id) const {
			return id <= max_id();
		}

		[[nodiscard]] size_t particle_count() const {
			return  particles.size();
		}

		[[nodiscard]] FetcherT get_fetcher(const size_t index) noexcept{
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			return FetcherT(particles[index]);
		}

	protected:
		void swap_particles (uint32_t idx1, uint32_t idx2) {
			auto id1 = particles[idx1].id, id2 = particles[idx2].id;
			std::swap(particles[idx1], particles[idx2]);
			std::swap(id_to_index_map[id1], id_to_index_map[id2]);
		}


		std::vector<Particle> particles = {};
		std::vector<uint32_t> id_to_index_map; // map id to index
	private:
		bool is_built = false;
	};
}