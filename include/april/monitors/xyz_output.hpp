#pragma once

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include "april/monitors/monitor.hpp"
#include "april/exec/kernel.hpp"


namespace april {

    class XYZOutput final : public monitor::Monitor {
    public:
        static constexpr auto fields =
            ParticleField::position |
            ParticleField::type;

        explicit XYZOutput(const Trigger & trigger,std::filesystem::path path = "trajectory.xyz"):
            Monitor(trigger),
            output_path(std::move(path))
        {}


        void initialize() {
            const std::filesystem::path parent = output_path.parent_path();
            if (!parent.empty()) {
                std::filesystem::create_directories(parent);
            }

            output.open(output_path, std::ios::out | std::ios::trunc);
            if (!output) {
                throw std::runtime_error("Failed to create XYZ trajectory: " + output_path.string());
            }

            output << std::setprecision(
                std::numeric_limits<double>::max_digits10
            );
        }


        template<class S>
        void record(const core::SystemContext<S> & sys) {
            if (!output.is_open()) {
                throw std::runtime_error(
                    "Cannot write XYZ trajectory because the output file is not open"
                );
            }

            output << sys.size() << '\n';
            output << "step=" << sys.step() << " time=" << sys.time() << '\n';

            sys.for_each_particle_view(scalar_kernel<fields>(
                [&](const auto & p) {
                    output
                        << 'T' << p.type << ' '
                        << p.position.x << ' '
                        << p.position.y << ' '
                        << p.position.z << '\n';
                }
            ));

            if (!output) {
                throw std::runtime_error("Failed to write XYZ trajectory: " + output_path.string());
            }
        }


        void finalize() {
            if (!output.is_open()) return;

            output.flush();
            if (!output) {
                throw std::runtime_error("Failed to flush XYZ trajectory: " + output_path.string());
            }

            output.close();
        }


    private:
        std::filesystem::path output_path;
        std::ofstream output;
    };

} // namespace april
