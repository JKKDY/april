#pragma once

#include "april/containers/container.h"

namespace april::cont::impl {
	template<typename Config, typename Env>
	class ContiguousContainer : public Container< Config, Env> {
	public:
		using Base = Container<Config, Env>;
		using typename Base::Particle;
		using typename Base::ParticleID;

		using Base::Base;

		void build_storage(const std::vector<Particle>& particles) {
			this->particles = std::vector(particles);
			for (size_t i = 0; i < particles.size(); i++) {
				indices.push_back(i);
			}
			this->is_built = true;
		}

		[[nodiscard]] Particle& get_particle_by_id(ParticleID id) noexcept{
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			AP_ASSERT(id >= id_start() && id <= id_end(), "invalid id. got " + std::to_string(id));
			return particles[indices[static_cast<size_t>(id)]];
		}

		[[nodiscard]] ParticleID id_start() const {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			return 0;
		}

		[[nodiscard]] ParticleID id_end() const {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			return particles.size() - 1; // note: will underflow if particles.size == 0
		}

		// index for non-stable iteration
		[[nodiscard]] Particle& get_particle_by_index(size_t index) noexcept {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			AP_ASSERT(index < particles.size(), "index must be < #particles");
			return particles[index];
		}

		[[nodiscard]] size_t index_start() const {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			return 0;
		}

		[[nodiscard]] size_t index_end() const {
			AP_ASSERT(is_built, "storage was not built. build_storage must be called");
			return particles.size() - 1; // note: will underflow if particles.size == 0
		}

		[[nodiscard]] size_t particle_count() const {
			return particles.size();
		}

	protected:
		void swap_particles (size_t idx1, size_t idx2) {
			std::swap(particles[idx1], particles[idx2]);
			std::swap(indices[idx1], indices[idx2]);
		}

		std::vector<Particle> particles;
		std::vector<size_t> indices;
	private:
		bool is_built = false;
	};
}