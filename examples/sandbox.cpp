#include <april/april.hpp>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

int main() {
    const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/sandbpox";
    remove_all(dir_path);   // delete the directory and all contents
    create_directory(dir_path); // recreate the empty directory

    auto liquid = ParticleCuboid()
        .at({1.5, 2, 0})
        .velocity({0, 0, 0})
        .count({50, 50, 1})
        .mass(1.0)
        .spacing(1.2)
        .type(0)
        .thermal([](vec3 /*pos*/) {
                constexpr double avg_vel = 1.0;
                 constexpr unsigned dim = 2;
                return math::maxwell_boltzmann_velocity(avg_vel, dim);
            });

    // auto drop = ParticleSphere()
    //     .at({20,20,0})
    //     .radius_xyz({10, 10, 0})
    //     .mass(1.0)
    //     .spacing(1)
    //     .type(1)
    //     .thermal([](vec3 /*pos*/) {
    //             constexpr double avg_vel = 1.0;
    //             return math::maxwell_boltzmann_velocity(avg_vel);
    //         });

    auto thermostat = VelocityScalingThermostat(0.5,
        temperature_not_set,
        temperature_not_set,
        Trigger::every(1000));

    auto gravity = UniformField({0, -12.44, 0});

    auto env = Environment(
        boundaries<ReflectiveBoundary>,
        forces<LennardJones>,
        controllers<VelocityScalingThermostat>,
        fields<UniformField>);

    env.with_extent(160,180, 0)
        .with_force(LennardJones(3, 1.2), to_type(0))
        // .with_force(LennardJones(3, 1), to_type(1))
        .with_particles(liquid)
        // .with_particles(drop)
        .with_boundaries(ReflectiveBoundary(), all_faces)
        .with_controller(thermostat)
        .with_field(gravity);

    auto container = LinkedCells<Layout::AoSoA<>>();
    auto system = build_system(env, container);

    auto integrator = VelocityVerlet(system, monitors<Benchmark, ProgressBar, BinaryOutput>)
        .with_monitor(Benchmark())
        .with_monitor(BinaryOutput(Trigger::every(100), dir_path.string()))
        .with_monitor(ProgressBar(Trigger::every(100)))
        .with_dt(0.0002)
        .for_duration(50)
        .run();
}
