#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstring>



# include "april/april.h"
using namespace april;
namespace fs = std::filesystem;


// helper to read raw bytes
template<typename T> static T read_binary(std::ifstream& in) {
	T val;
	in.read(reinterpret_cast<char*>(&val), sizeof(val));
	return val;
}

// helper to create dummy particle
static ext::Particle make_particle(env::impl::ParticleType type, env::impl::ParticleID id, vec3 pos={0,0,0}, ParticleState state= ParticleState::ALIVE) {
	return {
		/* id          */ id,
		/* position    */ pos,
		/* velocity    */ vec3{0,0,0},
		/* mass        */ 1.0,
		/* type        */ type,
		/* state       */ state
		};
}


class BinaryOutputTest : public testing::Test {
protected:
	fs::path dir;
	std::string base{"bin_test"};

	void SetUp() override {
		dir = fs::path(".") / "binary_output_test";
		fs::remove_all(dir);
		fs::create_directories(dir);
	}
	void TearDown() override {
		fs::remove_all(dir);
	}
};


// TEST 1: Header only, zero particles
TEST_F(BinaryOutputTest, EmptyFileContainsOnlyHeader) {
	std::vector<ParticleView> empty;
	BinaryOutput out(1, dir.string(), base);

	out.record(0, 0, empty);

	auto path = dir / (base + "_00000.bin");
	ASSERT_TRUE(fs::exists(path));

	std::ifstream in{path, std::ios::binary};
	ASSERT_TRUE(in.good());

	// Magic
	char magic[4]; in.read(magic,4);
	constexpr char expected_magic[4] = { 'P', 'A', 'R', 'T' };
	EXPECT_EQ(std::memcmp(magic, expected_magic, 4), 0);

	// Version
	EXPECT_EQ(read_binary<uint32_t>(in), 1u);

	// Step
	EXPECT_EQ(read_binary<uint64_t>(in), 0u);

	// Count = 0
	EXPECT_EQ(read_binary<uint64_t>(in), 0u);

	// Flags
	EXPECT_EQ(read_binary<uint32_t>(in), 0u);

	// No further data
	EXPECT_TRUE(in.peek() == EOF);
}


// TEST 2: Single particle record
TEST_F(BinaryOutputTest, SingleParticle) {
	auto p = make_particle(5, 2, vec3{1,2,3}, ParticleState::ALIVE);
	std::vector v {ParticleView(p)};
	BinaryOutput out(1, dir.string(), base);

	out.record(1, 0, v);

	auto path = dir / (base + "_00001.bin");
	std::ifstream in{path, std::ios::binary};
	// skip header fields (13 bytes after magic+version+step+count+flags)
	in.seekg(4 + 4 + 8 + 8 + 4);

	// Read 3 floats
	auto fx = read_binary<float>(in);
	auto fy = read_binary<float>(in);
	auto fz = read_binary<float>(in);
	EXPECT_FLOAT_EQ(fx, 1.0f);
	EXPECT_FLOAT_EQ(fy, 2.0f);
	EXPECT_FLOAT_EQ(fz, 3.0f);

	// type,id,state
	EXPECT_EQ(read_binary<uint32_t>(in), 5u);
	EXPECT_EQ(read_binary<uint32_t>(in), 2u);
	EXPECT_EQ(read_binary<uint8_t>(in), static_cast<uint8_t>(ParticleState::ALIVE));
}


// TEST 3: Multiple particle records
TEST_F(BinaryOutputTest, MultipleParticles) {
	auto p1 = make_particle(1, 0, vec3{0,0,0}, ParticleState::DEAD);
	auto p2 = make_particle(2, 1, vec3{4,5,6}, ParticleState::ALIVE);
	auto p3 = make_particle(3, 2, vec3{7,8,9}, ParticleState::PASSIVE);
	std::vector v {ParticleView(p1),ParticleView(p2),ParticleView(p3)};
	BinaryOutput out(1, dir.string(), base);

	out.record(2, 0, v);

	auto path = dir / (base + "_00002.bin");
	std::ifstream in{path, std::ios::binary};
	// skip header fields (13 bytes after magic+version+step+count+flags)
	in.seekg(4 + 4 + 8 + 8 + 4);

	for (auto & i : v) {
		auto x = read_binary<float>(in);
		auto y = read_binary<float>(in);
		auto z = read_binary<float>(in);
		EXPECT_FLOAT_EQ(x, static_cast<float>(i.position.x));
		EXPECT_FLOAT_EQ(y, static_cast<float>(i.position.y));
		EXPECT_FLOAT_EQ(z, static_cast<float>(i.position.z));

		EXPECT_EQ(read_binary<uint32_t>(in), i.type);
		EXPECT_EQ(read_binary<uint32_t>(in), i.id);
		EXPECT_EQ(read_binary<uint8_t>(in), static_cast<uint8_t>(i.state));
	}
}
