#include <april/april.hpp>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;




int main() {
    const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/sandbox";
    remove_all(dir_path);   // delete the directory and all contents
    create_directory(dir_path); // recreate the empty directory


    struct ChargeAttribute {
        double charge;

        struct VectorLayout {
            using ScalarType = double;
            simd::Packed<ScalarType> charge;
        };
    };


    int DEFAULT = 0;

    // 1) Define particles and interactions
    auto sun = Particle().at(0, 0, 0).with_mass(1.0).as_type(DEFAULT).with_data(ChargeAttribute{1});
    auto planet = Particle().at(1, 0, 0).with_velocity(0, 1, 0).with_mass(1e-3).as_type(DEFAULT).with_data(ChargeAttribute{1});
    auto moon = Particle().at(1.05, 0, 0).with_velocity(0, 1.2, 0).with_mass(1e-6).as_type(DEFAULT).with_data(ChargeAttribute{1});

    // Declare which component types may be used
    auto env = Environment(forces<Coulomb>, boundaries<OpenBoundary>, particle_attributes<ChargeAttribute>)
        .with_particles({sun, planet, moon})
        .with_force(Coulomb(), to_type(DEFAULT));

    constexpr auto algo = DirectSum<Layout::SoA>();
    auto system = build_system(env, algo);

    auto integrator = Yoshida4(system, monitors<BinaryOutput, ProgressBar, Benchmark, TerminalOutput>)
        .with_monitor(BinaryOutput(Trigger::every(50), dir_path.string()))
        .with_monitor(ProgressBar(Trigger::every(50)))
        .with_monitor(Benchmark())
        .run_for_duration(0.014, 1000);
}

