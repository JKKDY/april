#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <type_traits>
#include <utility>

#include "april/april.hpp"
#include "april/monitors/vtp_output.hpp"

using namespace april;
using namespace april::core;
using testing::HasSubstr;

namespace fs = std::filesystem;


// helper to read a text file
static std::string read_file(const fs::path & path) {
	std::ifstream in(path);
	return {
		std::istreambuf_iterator<char>(in),
		std::istreambuf_iterator<char>()
	};
}


// helper to count occurrences in a string
static size_t count_occurrences(
	const std::string & text,
	const std::string & value
) {
	size_t count = 0;
	size_t pos = text.find(value);

	while (pos != std::string::npos) {
		count++;
		pos = text.find(value, pos + value.size());
	}

	return count;
}


// helper to create dummy particle
using ParticleRec = particle::ParticleRecord<NoParticleAttributes>;

static ParticleRec make_particle_rec(
	const ParticleType type,
	const ParticleID id,
	const vec3 & position = {0,0,0},
	const vec3 & velocity = {0,0,0},
	const vec3 & force = {0,0,0}
) {
	ParticleRec rec;
	rec.id = id;
	rec.position = position;
	rec.velocity = velocity;
	rec.force = force;
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
		particle::internal::ScalarParticleRef<
			M,
			ParticleField::none,
			user_data_t
		>;

	explicit DummySystem(
		const size_t step,
		const double time,
		std::vector<ParticleRec> particle_data
	)
		: particles(std::move(particle_data))
		, step_(step)
		, time_(time)
	{}

	std::vector<ParticleRec> particles;
	size_t step_ = 0;
	double time_ = 0.0;
	Box sim_box = {{0,0,0}, {1,1,1}};

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
		Kernel&& func,
		ParticleState = ParticleState::ALL
	) const {
		using K = std::remove_cvref_t<Kernel>;

		for (size_t i = 0; i < size(); i++) {
			const auto& p = view<K::Read>(i);
			func(p);
		}
	}

	template<ParticleField M>
	[[nodiscard]] auto view(const size_t index) const noexcept {
		const ParticleRec& record = particles.at(index);

		auto get_field = [&]<ParticleField F>() {
			if constexpr (F == ParticleField::position)
				return &record.position;
			else if constexpr (F == ParticleField::velocity)
				return &record.velocity;
			else if constexpr (F == ParticleField::force)
				return &record.force;
			else if constexpr (F == ParticleField::old_position)
				return &record.old_position;
			else if constexpr (F == ParticleField::mass)
				return &record.mass;
			else if constexpr (F == ParticleField::state)
				return &record.state;
			else if constexpr (F == ParticleField::type)
				return &record.type;
			else if constexpr (F == ParticleField::id)
				return &record.id;
			else if constexpr (F == ParticleField::attributes)
				return &record.attributes;
		};

		auto src = particle::internal::make_particle_source<
			M,
			ParticleField::none
		>(get_field);

		return ParticleView<M>(src);
	}

	const SystemContext<DummySystem> ctx = SystemContext(*this);
};


class VTPOutputTest : public testing::Test {
protected:
	fs::path dir;
	std::string base{"particles"};

	void SetUp() override {
		dir = fs::path(".") / "vtp_output_test";
		fs::remove_all(dir);
	}

	void TearDown() override {
		fs::remove_all(dir);
	}
};


// TEST 1: initialize creates the requested output directory
TEST_F(VTPOutputTest, InitializeCreatesOutputDirectory) {
	const fs::path nested_dir = dir / "nested" / "output";

	VTPOutput<> out(
		Trigger::always(),
		nested_dir,
		base
	);

	out.initialize();

	EXPECT_TRUE(fs::exists(nested_dir));
	EXPECT_TRUE(fs::is_directory(nested_dir));
}


// TEST 2: finalize creates a valid empty collection
TEST_F(VTPOutputTest, FinalizeWithoutFramesCreatesEmptyCollection) {
	VTPOutput<> out(
		Trigger::always(),
		dir,
		base
	);

	out.initialize();
	out.finalize();

	const fs::path path = dir / (base + ".pvd");
	ASSERT_TRUE(fs::exists(path));

	const std::string xml = read_file(path);

	EXPECT_THAT(xml, HasSubstr("<?xml version=\"1.0\"?>"));
	EXPECT_THAT(xml, HasSubstr("<VTKFile type=\"Collection\""));
	EXPECT_THAT(xml, HasSubstr("<Collection>"));
	EXPECT_THAT(xml, HasSubstr("</Collection>"));
	EXPECT_THAT(xml, HasSubstr("</VTKFile>"));

	EXPECT_EQ(
		count_occurrences(xml, "<DataSet"),
		0u
	);
}


// TEST 3: empty systems still produce a valid VTP frame
TEST_F(VTPOutputTest, EmptySystemProducesValidFrame) {
	VTPOutput<> out(
		Trigger::always(),
		dir,
		base
	);

	DummySystem sys(0, 0.0, {});

	out.initialize();
	out.record(sys.ctx);
	out.finalize();

	const fs::path frame_path =
		dir /
		(base + "_00000.vtp");

	ASSERT_TRUE(fs::exists(frame_path));

	const std::string xml = read_file(frame_path);

	EXPECT_THAT(xml, HasSubstr("<VTKFile type=\"PolyData\""));
	EXPECT_THAT(xml, HasSubstr("<PolyData>"));

	EXPECT_THAT(
		xml,
		HasSubstr(
			"<Piece NumberOfPoints=\"0\" "
			"NumberOfVerts=\"0\" "
			"NumberOfLines=\"0\" "
			"NumberOfStrips=\"0\" "
			"NumberOfPolys=\"0\">"
		)
	);

	EXPECT_THAT(xml, HasSubstr("<PointData>"));
	EXPECT_THAT(xml, HasSubstr("<Points>"));
	EXPECT_THAT(xml, HasSubstr("<Verts>"));
	EXPECT_THAT(xml, HasSubstr("Name=\"connectivity\""));
	EXPECT_THAT(xml, HasSubstr("Name=\"offsets\""));
}


// TEST 4: default output writes position, velocity, force, type and ID
TEST_F(VTPOutputTest, WritesDefaultParticleFields) {
	std::vector particles {
		make_particle_rec(
			2,
			10,
			vec3{13,14,15},
			vec3{1,2,3},
			vec3{7,8,9}
		),
		make_particle_rec(
			5,
			20,
			vec3{16,17,18},
			vec3{4,5,6},
			vec3{10,11,12}
		)
	};

	VTPOutput<> out(
		Trigger::always(),
		dir,
		base
	);

	DummySystem sys(4, 0.125, std::move(particles));

	out.initialize();
	out.record(sys.ctx);
	out.finalize();

	const fs::path frame_path =
		dir /
		(base + "_00004.vtp");

	ASSERT_TRUE(fs::exists(frame_path));

	const std::string xml = read_file(frame_path);

	EXPECT_THAT(
		xml,
		HasSubstr(
			"<Piece NumberOfPoints=\"2\" "
			"NumberOfVerts=\"2\""
		)
	);

	EXPECT_THAT(
		xml,
		HasSubstr(
			"<DataArray type=\"Float64\" "
			"Name=\"velocity\" "
			"NumberOfComponents=\"3\" "
			"format=\"ascii\">"
		)
	);

	EXPECT_THAT(
		xml,
		HasSubstr(
			"<DataArray type=\"Float64\" "
			"Name=\"force\" "
			"NumberOfComponents=\"3\" "
			"format=\"ascii\">"
		)
	);

	EXPECT_THAT(
		xml,
		HasSubstr(
			"<DataArray type=\"UInt32\" "
			"Name=\"type\" "
			"format=\"ascii\">"
		)
	);

	EXPECT_THAT(
		xml,
		HasSubstr(
			"<DataArray type=\"UInt32\" "
			"Name=\"id\" "
			"format=\"ascii\">"
		)
	);

	EXPECT_THAT(xml, HasSubstr("1 2 3\n4 5 6\n"));
	EXPECT_THAT(xml, HasSubstr("7 8 9\n10 11 12\n"));
	EXPECT_THAT(xml, HasSubstr("13 14 15\n16 17 18\n"));
	EXPECT_THAT(xml, HasSubstr("2\n5\n"));
	EXPECT_THAT(xml, HasSubstr("10\n20\n"));
}


// TEST 5: compile-time field selection excludes disabled fields
TEST_F(VTPOutputTest, WritesOnlyConfiguredFields) {
	using CompactOutput = VTPOutput<
		ParticleField::position |
		ParticleField::type
	>;

	std::vector particles {
		make_particle_rec(
			3,
			10,
			vec3{1,2,3},
			vec3{4,5,6},
			vec3{7,8,9}
		)
	};

	CompactOutput out(
		Trigger::always(),
		dir,
		base
	);

	DummySystem sys(2, 0.25, std::move(particles));

	out.initialize();
	out.record(sys.ctx);
	out.finalize();

	const fs::path frame_path =
		dir /
		(base + "_00002.vtp");

	ASSERT_TRUE(fs::exists(frame_path));

	const std::string xml = read_file(frame_path);

	EXPECT_THAT(xml, HasSubstr("Name=\"type\""));

	EXPECT_EQ(
		xml.find("Name=\"velocity\""),
		std::string::npos
	);

	EXPECT_EQ(
		xml.find("Name=\"force\""),
		std::string::npos
	);

	EXPECT_EQ(
		xml.find("Name=\"id\""),
		std::string::npos
	);

	EXPECT_THAT(xml, HasSubstr("1 2 3\n"));
}


// TEST 6: each particle is represented by one vertex cell
TEST_F(VTPOutputTest, WritesVertexConnectivityAndOffsets) {
	std::vector particles {
		make_particle_rec(0, 0),
		make_particle_rec(0, 1),
		make_particle_rec(0, 2)
	};

	VTPOutput<> out(
		Trigger::always(),
		dir,
		base
	);

	DummySystem sys(3, 0.0, std::move(particles));

	out.initialize();
	out.record(sys.ctx);
	out.finalize();

	const fs::path frame_path =
		dir /
		(base + "_00003.vtp");

	ASSERT_TRUE(fs::exists(frame_path));

	const std::string xml = read_file(frame_path);

	const size_t connectivity_pos =
		xml.find("Name=\"connectivity\"");

	ASSERT_NE(connectivity_pos, std::string::npos);

	const size_t connectivity_end =
		xml.find("</DataArray>", connectivity_pos);

	ASSERT_NE(connectivity_end, std::string::npos);

	const std::string connectivity =
		xml.substr(
			connectivity_pos,
			connectivity_end - connectivity_pos
		);

	EXPECT_THAT(connectivity, HasSubstr("0\n1\n2\n"));


	const size_t offsets_pos =
		xml.find("Name=\"offsets\"");

	ASSERT_NE(offsets_pos, std::string::npos);

	const size_t offsets_end =
		xml.find("</DataArray>", offsets_pos);

	ASSERT_NE(offsets_end, std::string::npos);

	const std::string offsets =
		xml.substr(
			offsets_pos,
			offsets_end - offsets_pos
		);

	EXPECT_THAT(offsets, HasSubstr("1\n2\n3\n"));
}


// TEST 7: multiple frames are added to the PVD collection
TEST_F(VTPOutputTest, WritesMultipleFramesToCollection) {
	std::vector particles {
		make_particle_rec(1, 0, vec3{0,0,0})
	};

	VTPOutput<> out(
		Trigger::always(),
		dir,
		base
	);

	DummySystem sys(5, 0.125, std::move(particles));

	out.initialize();
	out.record(sys.ctx);

	sys.step_ = 10;
	sys.time_ = 0.25;
	sys.particles[0].position = vec3{1,2,3};

	out.record(sys.ctx);
	out.finalize();

	const fs::path first_frame =
		dir /
		(base + "_00005.vtp");

	const fs::path second_frame =
		dir /
		(base + "_00010.vtp");

	const fs::path collection =
		dir /
		(base + ".pvd");

	ASSERT_TRUE(fs::exists(first_frame));
	ASSERT_TRUE(fs::exists(second_frame));
	ASSERT_TRUE(fs::exists(collection));

	const std::string xml = read_file(collection);

	EXPECT_EQ(
		count_occurrences(xml, "<DataSet"),
		2u
	);

	EXPECT_THAT(
		xml,
		HasSubstr(
			"timestep=\"0.125\" "
			"group=\"\" "
			"part=\"0\" "
			"file=\"particles_00005.vtp\""
		)
	);

	EXPECT_THAT(
		xml,
		HasSubstr(
			"timestep=\"0.25\" "
			"group=\"\" "
			"part=\"0\" "
			"file=\"particles_00010.vtp\""
		)
	);
}


// TEST 8: frame files use the configured base name
TEST_F(VTPOutputTest, UsesConfiguredBaseName) {
	const std::string custom_base = "trajectory";

	VTPOutput<> out(
		Trigger::always(),
		dir,
		custom_base
	);

	DummySystem sys(
		12,
		0.5,
		{make_particle_rec(0, 0)}
	);

	out.initialize();
	out.record(sys.ctx);
	out.finalize();

	EXPECT_TRUE(
		fs::exists(
			dir /
			"trajectory_00012.vtp"
		)
	);

	EXPECT_TRUE(
		fs::exists(
			dir /
			"trajectory.pvd"
		)
	);
}