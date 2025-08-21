#pragma once


#include "april/algo/algorithm.h"

namespace april::algo::impl {
	class DirectSum;
}

namespace april::algo {

	class DirectSum {
		using impl = impl::DirectSum;
	};

	namespace impl {
		class DirectSum final : public Algorithm<algo::DirectSum> {
		public:
			using Algorithm::Algorithm;

			void build(const std::vector<Particle> & particles) override;
			void calculate_forces() override;
		private:
			std::vector<Particle> particles;
		};
	}
}