#pragma once

#include <gtest/gtest.h>
# include "april/april.h"

using namespace april;


class OrbitMonitor final : public Monitor {
public:
	OrbitMonitor(): Monitor(Trigger::always()) {}
	explicit OrbitMonitor(const double v, const double r): Monitor(Trigger::always()), v(v), r(r) {}

	static constexpr env::FieldMask fields = to_field_mask(env::Field::all);

	template<class S>
	void record(const SystemContext<S> & sys) const {
		EXPECT_EQ(sys.size(), 2u);
		const size_t id1 = sys.index_start();
		const size_t id2 = sys.index_end() - 1;
		const ParticleView p = sys.template get_particle_by_index<fields>(id1).mass < 1 ?
			sys.template get_particle_by_index<fields>(id1) : sys.template get_particle_by_index<fields>(id2);

		EXPECT_NEAR(p.velocity.norm(), v, 1e-3);
		EXPECT_NEAR(p.position.norm(), r, 1e-3);
	}

	double v{};
	double r{};
};
