#pragma once


#include "april/containers/container.h"

namespace april::cont::impl {
	class DirectSum;
}

namespace april::cont {

	struct DirectSum {
		using impl = impl::DirectSum;
	};

	namespace impl {
		class DirectSum final : public Container<cont::DirectSum> {
		public:
			using Container::Container;

			void build(const std::vector<Particle> & particles);
			void calculate_forces();

			[[nodiscard]] Particle & get_particle_by_id(ParticleID id);
			[[nodiscard]] ParticleID id_start() const;
			[[nodiscard]] ParticleID id_end() const;

			[[nodiscard]] Particle & get_particle_by_index(size_t index) noexcept;
			[[nodiscard]] size_t index_start() const;
			[[nodiscard]] size_t index_end() const;

			[[nodiscard]] size_t particle_count() const;
		private:
			std::vector<Particle> particles;
		};
	}
}