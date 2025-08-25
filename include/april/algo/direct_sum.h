#pragma once


#include "april/algo/algorithm.h"

namespace april::algo::impl {
	class DirectSum;
}

namespace april::algo {

	struct DirectSum {
		using impl = impl::DirectSum;
	};

	namespace impl {
		class DirectSum final : public Algorithm<algo::DirectSum> {
		public:
			using Algorithm::Algorithm;

			void build(const std::vector<Particle> & particles) override;
			void calculate_forces() override;

			Particle & get_particle_by_id(ParticleID id) override;
			ParticleID id_start() override;
			ParticleID id_end() override;

			Particle & get_particle_by_index(size_t index) noexcept override;
			size_t index_start() override;
			size_t index_end() override;

			size_t particle_count() override;
		private:
			std::vector<Particle> particles;
		};
	}
}