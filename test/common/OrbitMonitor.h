#pragma once

#include <gtest/gtest.h>
# include "april/april.h"

using namespace april;


class OrbitMonitor final : public Monitor {
public:
	OrbitMonitor(): Monitor(1) {}
	explicit OrbitMonitor(const double v, const double r): Monitor(1), v(v), r(r) {}

	void record(const SimulationContext & sys) const {
		EXPECT_EQ(sys.size(), 2u);
		const size_t id1 = sys.index_start();
		const size_t id2 = sys.index_end() - 1;
		const ParticleView p = sys.get_particle_by_index(id1).mass < 1 ?
			sys.get_particle_by_index(id1) : sys.get_particle_by_index(id2);

		EXPECT_NEAR(p.velocity.norm(), v, 1e-3);
		EXPECT_NEAR(p.position.norm(), r, 1e-3);
	}

	double v{};
	double r{};
};
