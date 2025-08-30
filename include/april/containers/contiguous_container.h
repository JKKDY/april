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

		[[nodiscard]] Particle& get_particle_by_id(ParticleID) {
			throw std::runtime_error("Not implemented yet");
		}

		[[nodiscard]] ParticleID id_start() const {
			return 0;
		}

		[[nodiscard]] ParticleID id_end() const {
			return particles.size() - 1;
		}

		[[nodiscard]] Particle& get_particle_by_index(size_t index) noexcept {
			AP_ASSERT(index < particles.size(), "index must be < #particles");
			return particles[index];
		}

		[[nodiscard]] size_t index_start() const {
			return 0;
		}

		[[nodiscard]] size_t index_end() const {
			return particles.size() - 1;
		}

		[[nodiscard]] size_t particle_count() const {
			return particles.size();
		}

	protected:
		std::vector<Particle> particles;

	};
}