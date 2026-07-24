#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "april/april.hpp"
#include "april/monitors/xyz_output.hpp"


using namespace april;
using namespace april::core;

namespace fs = std::filesystem;


using ParticleRec = particle::ParticleRecord<NoParticleAttributes>;

static ParticleRec make_particle_rec(
	const ParticleType type,
	const ParticleID id,
	const vec3 & position = {0, 0, 0}
) {
	ParticleRec rec;
	rec.id = id;
	rec.position = position;
	rec.velocity = vec3{0, 0, 0};
	rec.mass = 1.0;
	rec.type = type;
	rec.state = ParticleState::ALIVE;
	return rec;
}


class DummySystem {
public:
	using user_data_t = NoParticleAttributes;

	template<ParticleField M>
	using ParticleRef =
		particle::internal::ScalarParticleRef<M, M, user_data_t>;

	template<ParticleField M>
	using ParticleView =
		particle::internal::ScalarParticleRef<M, ParticleField::none, user_data_t>;


	explicit DummySystem(
		const size_t step,
		const double time,
		std::vector<ParticleRec> particle_data
	)
		:
		particles(std::move(particle_data)),
		step_(step),
		time_(time)
	{}


	std::vector<ParticleRec> particles;
	size_t step_ = 0;
	double time_ = 0.0;
	Box sim_box = {{0, 0, 0}, {1, 1, 1}};


	[[nodiscard]] size_t size(
		ParticleState = ParticleState::ALL
	) const noexcept {
		return particles.size();
	}

	[[nodiscard]] size_t step() const noexcept {
		return step_;
	}

	[[nodiscard]] double time() const noexcept {
		return time_;
	}

	[[nodiscard]] Box box() const noexcept {
		return sim_box;
	}


	template<
		ParallelPolicy P = ParallelPolicy::Serial,
		VectorPolicy V = VectorPolicy::Auto,
		typename Kernel
	>
	void for_each_particle_view(
		Kernel && func,
		ParticleState = ParticleState::ALL
	) const {
		using K = std::remove_cvref_t<Kernel>;

		for (size_t i = 0; i < size(); ++i) {
			const auto & p = view<K::Read>(i);
			func(p);
		}
	}


	template<ParticleField M>
	[[nodiscard]] auto view(const size_t index) const noexcept {
		const ParticleRec & record = particles.at(index);

		particle::internal::ParticleSource<
			M,
			ParticleField::none,
			user_data_t
		> src;

		if constexpr (
			particle::internal::has_field_v<M, ParticleField::position>
		) {
			src.position = &record.position;
		}

		if constexpr (
			particle::internal::has_field_v<M, ParticleField::velocity>
		) {
			src.velocity = &record.velocity;
		}

		if constexpr (
			particle::internal::has_field_v<M, ParticleField::force>
		) {
			src.force = &record.force;
		}

		if constexpr (
			particle::internal::has_field_v<M, ParticleField::old_position>
		) {
			src.old_position = &record.old_position;
		}

		if constexpr (
			particle::internal::has_field_v<M, ParticleField::mass>
		) {
			src.mass = &record.mass;
		}

		if constexpr (
			particle::internal::has_field_v<M, ParticleField::state>
		) {
			src.state = &record.state;
		}

		if constexpr (
			particle::internal::has_field_v<M, ParticleField::type>
		) {
			src.type = &record.type;
		}

		if constexpr (
			particle::internal::has_field_v<M, ParticleField::id>
		) {
			src.id = &record.id;
		}

		if constexpr (
			particle::internal::has_field_v<M, ParticleField::attributes>
		) {
			src.attributes = &record.attributes;
		}

		return ParticleView<M>(src);
	}


	const SystemContext<DummySystem> ctx = SystemContext(*this);
};


class XYZOutputTest : public testing::Test {
protected:
	fs::path dir;

	void SetUp() override {
		dir = fs::path(".") / "xyz_output_test";
		fs::remove_all(dir);
		fs::create_directories(dir);
	}

	void TearDown() override {
		fs::remove_all(dir);
	}
};


// TEST 1: Empty systems still produce a valid XYZ frame
TEST_F(XYZOutputTest, EmptySystemProducesValidFrame) {
	const fs::path path = dir / "nested" / "trajectory.xyz";

	XYZOutput out(Trigger::always(), path);
	DummySystem sys(0, 0.0, {});

	out.initialize();
	out.record(sys.ctx);
	out.finalize();

	ASSERT_TRUE(fs::exists(path));

	std::ifstream in(path);
	ASSERT_TRUE(in.good());

	std::string line;

	ASSERT_TRUE(std::getline(in, line));
	EXPECT_EQ(line, "0");

	ASSERT_TRUE(std::getline(in, line));
	EXPECT_EQ(line, "step=0 time=0");

	EXPECT_FALSE(std::getline(in, line));
}


// TEST 2: Particle types and positions are written correctly
TEST_F(XYZOutputTest, WritesParticleData) {
	const fs::path path = dir / "particles.xyz";

	std::vector particles {
		make_particle_rec(2, 0, vec3{1.0, 2.0, 3.0}),
		make_particle_rec(5, 1, vec3{-4.5, 6.25, 8.0})
	};

	XYZOutput out(Trigger::always(), path);
	DummySystem sys(4, 0.125, std::move(particles));

	out.initialize();
	out.record(sys.ctx);
	out.finalize();

	std::ifstream in(path);
	ASSERT_TRUE(in.good());

	size_t particle_count = 0;
	in >> particle_count;
	EXPECT_EQ(particle_count, 2u);

	std::string metadata;
	std::getline(in, metadata);
	std::getline(in, metadata);
	EXPECT_EQ(metadata, "step=4 time=0.125");

	std::string label;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	ASSERT_TRUE(in >> label >> x >> y >> z);
	EXPECT_EQ(label, "T2");
	EXPECT_DOUBLE_EQ(x, 1.0);
	EXPECT_DOUBLE_EQ(y, 2.0);
	EXPECT_DOUBLE_EQ(z, 3.0);

	ASSERT_TRUE(in >> label >> x >> y >> z);
	EXPECT_EQ(label, "T5");
	EXPECT_DOUBLE_EQ(x, -4.5);
	EXPECT_DOUBLE_EQ(y, 6.25);
	EXPECT_DOUBLE_EQ(z, 8.0);
}


// TEST 3: Multiple record calls append frames to the same trajectory
TEST_F(XYZOutputTest, WritesMultipleFrames) {
	const fs::path path = dir / "trajectory.xyz";

	std::vector particles {
		make_particle_rec(1, 0, vec3{0.0, 0.0, 0.0})
	};

	XYZOutput out(Trigger::always(), path);
	DummySystem sys(1, 0.1, std::move(particles));

	out.initialize();
	out.record(sys.ctx);

	sys.step_ = 2;
	sys.time_ = 0.2;
	sys.particles[0].position = vec3{1.0, 2.0, 3.0};

	out.record(sys.ctx);
	out.finalize();

	std::ifstream in(path);
	ASSERT_TRUE(in.good());

	size_t particle_count = 0;
	std::string metadata;
	std::string label;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	ASSERT_TRUE(in >> particle_count);
	EXPECT_EQ(particle_count, 1u);

	std::getline(in, metadata);
	std::getline(in, metadata);

	EXPECT_TRUE(metadata.starts_with("step=1 time="));
	EXPECT_DOUBLE_EQ(
		std::stod(metadata.substr(metadata.find("time=") + 5)),
		0.1
	);

	ASSERT_TRUE(in >> label >> x >> y >> z);
	EXPECT_EQ(label, "T1");
	EXPECT_DOUBLE_EQ(x, 0.0);
	EXPECT_DOUBLE_EQ(y, 0.0);
	EXPECT_DOUBLE_EQ(z, 0.0);

	ASSERT_TRUE(in >> particle_count);
	EXPECT_EQ(particle_count, 1u);

	std::getline(in, metadata);
	std::getline(in, metadata);

	EXPECT_TRUE(metadata.starts_with("step=2 time="));
	EXPECT_DOUBLE_EQ(
		std::stod(metadata.substr(metadata.find("time=") + 5)),
		0.2
	);

	ASSERT_TRUE(in >> label >> x >> y >> z);
	EXPECT_EQ(label, "T1");
	EXPECT_DOUBLE_EQ(x, 1.0);
	EXPECT_DOUBLE_EQ(y, 2.0);
	EXPECT_DOUBLE_EQ(z, 3.0);

	in >> std::ws;
	EXPECT_TRUE(in.eof());
}