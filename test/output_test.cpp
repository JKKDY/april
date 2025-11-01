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
static env::internal::Particle make_particle(env::internal::ParticleType type, env::internal::ParticleID id, vec3 pos={0,0,0}, ParticleState state= ParticleState::ALIVE) {
	return {
		/* id          */ id,
		/* position    */ pos,
		/* velocity    */ vec3{0,0,0},
		/* mass        */ 1.0,
		/* type        */ type,
		/* state       */ state
		};
}



class DummyContext final : public SystemContext {
public:
    using ParticleView = env::ParticleView;
    using ParticleRef  = env::ParticleRef;
    using ParticleID   = env::internal::ParticleID;

    explicit DummyContext(const size_t step = 0, const double time = 0.0,
                          std::vector<ParticleView> particles = {})
        : step_(step), time_(time), particles_(std::move(particles)) {}
	env::Domain dummy_domain = {{0, 0, 0}, {1, 1, 1}};
    // ---- Core information ----
    [[nodiscard]] env::Domain domain() const noexcept override {
        return dummy_domain;
    }

	[[nodiscard]] env::Box box() const noexcept override {
    	return env::Box::from_domain(dummy_domain);
    }

    [[nodiscard]] double time() const noexcept override { return time_; }
    [[nodiscard]] size_t step() const noexcept override { return step_; }

    // ---- Particle access / modification ----
    [[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Box&) const override { return {}; }
    [[nodiscard]] std::vector<size_t> collect_indices_in_region(const env::Domain&) const override { return {}; }
    void register_particle_movement(ParticleID) override {}
    void register_all_particle_movements() override {}

    // ---- ID and index access ----

    [[nodiscard]] ParticleID id_start() const noexcept override { return 0; }
    [[nodiscard]] ParticleID id_end() const noexcept override { return 0; }


    [[nodiscard]] size_t index_start() const noexcept override { return 0; }
    [[nodiscard]] size_t index_end() const noexcept override { return particles_.size(); }

    // helper for test data
    [[nodiscard]] const std::vector<ParticleView>& particles() const noexcept { return particles_; }

	[[nodiscard]] ParticleView get_particle_by_index(const size_t idx) const noexcept override {
    	return particles_[idx];
    }

	[[nodiscard]] size_t size() const noexcept override {return particles_.size();}
	[[nodiscard]] size_t size(ParticleState ) const noexcept override {return particles_.size();}

private:
	[[nodiscard]] ParticleRef get_particle_by_id(ParticleID) noexcept override { std::unreachable(); }
	[[nodiscard]] ParticleView get_particle_by_id(ParticleID) const noexcept override { std::unreachable(); }
	[[nodiscard]] ParticleRef get_particle_by_index(size_t) noexcept override { std::unreachable();}

    size_t step_;
    double time_;
    std::vector<ParticleView> particles_;
};



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
	BinaryOutput out(Trigger::always(), dir.string(), base);

	DummyContext ctx(0, 0.0, empty);
	out.record(ctx);

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
	BinaryOutput out(Trigger::always(), dir.string(), base);

	DummyContext ctx(1, 0.0, v);
	out.record(ctx);

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
	BinaryOutput out(Trigger::always(), dir.string(), base);

	DummyContext ctx(2, 0.0, v);
	out.record(ctx);

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
