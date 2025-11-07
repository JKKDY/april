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
using ParticleRec =  env::internal::ParticleRecord<env::NoUserData>;
static ParticleRec make_particle_rec(ParticleType type, ParticleID id, vec3 pos={0,0,0}, ParticleState state= ParticleState::ALIVE) {
	ParticleRec rec;
	rec.id = id,
	rec.position = pos,
	rec.velocity = vec3{0,0,0},
	rec.mass = 1.0,
	rec.type = type,
	rec.state = state;
	return rec;
}

using namespace april;
using namespace april::env;

class DummySystem {
public:
	using user_data_t = NoUserData;
	template<FieldMask M> using ParticleRef         = ParticleRef<M, user_data_t>;
	template<FieldMask M> using ParticleView        = ParticleView<M, user_data_t>;
	template<FieldMask M> using RestrictedParticleRef = RestrictedParticleRef<M, user_data_t>;

	explicit DummySystem(
		size_t step,
		double time,
		std::vector<ParticleRec> particle_data
	) : particles(std::move(particle_data)), step_(step), time_(time)
	{}

	// mock data
	std::vector<ParticleRec> particles;
	size_t step_ = 0;
	double time_ = 0.0;
	Box sim_box = {{0,0,0}, {1,1,1}};

	// Implement minimal API for read-only tests
	[[nodiscard]] size_t index_start() const noexcept { return 0; }
	[[nodiscard]] size_t index_end() const noexcept { return particles.size(); }
	[[nodiscard]] size_t size(ParticleState = ParticleState::ALL) const noexcept { return particles.size(); }
	[[nodiscard]] size_t step() const noexcept { return step_; }
	[[nodiscard]] double time() const noexcept { return time_; }
	[[nodiscard]] Box box() const noexcept { return sim_box; }

	template<FieldMask M>
	[[nodiscard]] ParticleView<M> get_particle_by_index(const size_t index) const noexcept {
		const ParticleRec& record = particles.at(index);
		internal::ConstParticleRecordFetcher fetcher(record);
		return ParticleView<M>(fetcher);
	}

	const SystemContext<DummySystem> ctx = SystemContext<DummySystem>(*this);
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
	std::vector<ParticleRec> empty;
	BinaryOutput out(Trigger::always(), dir.string(), base);

	DummySystem sys(0, 0.0, empty);
	out.record(sys.ctx);

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
	auto p = make_particle_rec(5, 2, vec3{1,2,3}, ParticleState::ALIVE);
	std::vector v {ParticleRec(p)};
	BinaryOutput out(Trigger::always(), dir.string(), base);

	DummySystem sys(1, 0.0, v);
	out.record(sys.ctx);

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
	auto p1 = make_particle_rec(1, 0, vec3{0,0,0}, ParticleState::DEAD);
	auto p2 = make_particle_rec(2, 1, vec3{4,5,6}, ParticleState::ALIVE);
	auto p3 = make_particle_rec(3, 2, vec3{7,8,9}, ParticleState::PASSIVE);
	std::vector v {ParticleRec(p1),ParticleRec(p2),ParticleRec(p3)};
	BinaryOutput out(Trigger::always(), dir.string(), base);

	DummySystem sys(2, 0.0, v);
	out.record(sys.ctx);

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
