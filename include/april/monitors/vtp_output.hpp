#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include "april/monitors/monitor.hpp"
#include "april/utility/xml.hpp"


namespace april {

	template<
		ParticleField Fields =
			ParticleField::position |
			ParticleField::velocity |
			ParticleField::force |
			ParticleField::type |
			ParticleField::id
	>
	class VTPOutput final : public monitor::Monitor {
	public:
		static_assert(
			particle::internal::has_field_v<
				Fields,
				ParticleField::position
			>,
			"VTPOutput requires ParticleField::position"
		);

		static constexpr ParticleField fields = Fields;


		explicit VTPOutput(
			const Trigger & trigger,
			std::filesystem::path directory = ".",
			std::string base_name = "particles"
		)
			:
			Monitor(trigger),
			output_directory(std::move(directory)),
			base_name(std::move(base_name))
		{}


		void initialize() {
			std::error_code error;

			std::filesystem::create_directories(
				output_directory,
				error
			);

			if (error) {
				throw std::runtime_error(
					"Failed to create VTP output directory: " +
					output_directory.string()
				);
			}

			frames.clear();

			collection_path =
				output_directory /
				(base_name + ".pvd");
		}


		template<class S>
		void record(const core::SystemContext<S> & sys) {
			const std::size_t particle_count = sys.size();

			std::ostringstream file_name;
			file_name
				<< base_name
				<< '_'
				<< std::setfill('0')
				<< std::setw(5)
				<< sys.step()
				<< ".vtp";

			const std::filesystem::path relative_path = file_name.str();
			const std::filesystem::path frame_path =
				output_directory /
				relative_path;

			std::ofstream output(
				frame_path,
				std::ios::out |
				std::ios::trunc
			);

			if (!output) {
				throw std::runtime_error(
					"Failed to create VTP output file: " +
					frame_path.string()
				);
			}

			output << std::setprecision(
				std::numeric_limits<double>::max_digits10
			);

			utility::XMLWriter xml(output);

			xml.declaration();

			xml.open(
				"VTKFile",
				utility::attribute("type", "PolyData"),
				utility::attribute("version", "0.1"),
				utility::attribute("byte_order", "LittleEndian")
			);

			xml.open("PolyData");

			xml.open(
				"Piece",
				utility::attribute("NumberOfPoints", particle_count),
				utility::attribute("NumberOfVerts", particle_count),
				utility::attribute("NumberOfLines", 0),
				utility::attribute("NumberOfStrips", 0),
				utility::attribute("NumberOfPolys", 0)
			);


			// Particle fields
			xml.open("PointData");

			if constexpr (
				particle::internal::has_field_v<
					Fields,
					ParticleField::velocity
				>
			) {
				xml.open(
					"DataArray",
					utility::attribute("type", "Float64"),
					utility::attribute("Name", "velocity"),
					utility::attribute("NumberOfComponents", 3),
					utility::attribute("format", "ascii")
				);

				sys.for_each_particle_view(scalar_kernel<fields>(
					[&](const auto & p) {
						output
							<< p.velocity.x << ' '
							<< p.velocity.y << ' '
							<< p.velocity.z << '\n';
					}
				));

				xml.close("DataArray");
			}

			if constexpr (
				particle::internal::has_field_v<
					Fields,
					ParticleField::force
				>
			) {
				xml.open(
					"DataArray",
					utility::attribute("type", "Float64"),
					utility::attribute("Name", "force"),
					utility::attribute("NumberOfComponents", 3),
					utility::attribute("format", "ascii")
				);

				sys.for_each_particle_view(scalar_kernel<fields>(
					[&](const auto & p) {
						output
							<< p.force.x << ' '
							<< p.force.y << ' '
							<< p.force.z << '\n';
					}
				));

				xml.close("DataArray");
			}

			if constexpr (
				particle::internal::has_field_v<
					Fields,
					ParticleField::type
				>
			) {
				xml.open(
					"DataArray",
					utility::attribute("type", "UInt32"),
					utility::attribute("Name", "type"),
					utility::attribute("format", "ascii")
				);

				sys.for_each_particle_view(scalar_kernel<fields>(
					[&](const auto & p) {
						output
							<< static_cast<std::uint32_t>(p.type)
							<< '\n';
					}
				));

				xml.close("DataArray");
			}

			if constexpr (
				particle::internal::has_field_v<
					Fields,
					ParticleField::id
				>
			) {
				xml.open(
					"DataArray",
					utility::attribute("type", "UInt32"),
					utility::attribute("Name", "id"),
					utility::attribute("format", "ascii")
				);

				sys.for_each_particle_view(scalar_kernel<fields>(
					[&](const auto & p) {
						output
							<< static_cast<std::uint32_t>(p.id)
							<< '\n';
					}
				));

				xml.close("DataArray");
			}

			xml.close("PointData");


			// Particle positions
			xml.open("Points");

			xml.open(
				"DataArray",
				utility::attribute("type", "Float64"),
				utility::attribute("NumberOfComponents", 3),
				utility::attribute("format", "ascii")
			);

			sys.for_each_particle_view(scalar_kernel<fields>(
				[&](const auto & p) {
					output
						<< p.position.x << ' '
						<< p.position.y << ' '
						<< p.position.z << '\n';
				}
			));

			xml.close("DataArray");
			xml.close("Points");


			// One vertex cell per particle
			xml.open("Verts");

			xml.open(
				"DataArray",
				utility::attribute("type", "Int64"),
				utility::attribute("Name", "connectivity"),
				utility::attribute("format", "ascii")
			);

			for (std::size_t i = 0; i < particle_count; i++) {
				output
					<< static_cast<std::int64_t>(i)
					<< '\n';
			}

			xml.close("DataArray");

			xml.open(
				"DataArray",
				utility::attribute("type", "Int64"),
				utility::attribute("Name", "offsets"),
				utility::attribute("format", "ascii")
			);

			for (std::size_t i = 0; i < particle_count; i++) {
				output
					<< static_cast<std::int64_t>(i + 1)
					<< '\n';
			}

			xml.close("DataArray");
			xml.close("Verts");


			xml.close("Piece");
			xml.close("PolyData");
			xml.close("VTKFile");

			output.flush();

			if (!output) {
				throw std::runtime_error(
					"Failed to write VTP output file: " +
					frame_path.string()
				);
			}

			frames.push_back({
				sys.time(),
				relative_path
			});
		}


		void finalize() {
			std::ofstream output(
				collection_path,
				std::ios::out |
				std::ios::trunc
			);

			if (!output) {
				throw std::runtime_error(
					"Failed to create PVD collection file: " +
					collection_path.string()
				);
			}

			output << std::setprecision(std::numeric_limits<double>::max_digits10);

			utility::XMLWriter xml(output);

			xml.declaration();

			xml.open(
				"VTKFile",
				utility::attribute("type", "Collection"),
				utility::attribute("version", "0.1"),
				utility::attribute("byte_order", "LittleEndian")
			);

			xml.open("Collection");

			for (const auto & frame : frames) {
				xml.empty(
					"DataSet",
					utility::attribute("timestep", frame.time),
					utility::attribute("group", ""),
					utility::attribute("part", 0),
					utility::attribute("file", frame.file.generic_string())
				);
			}

			xml.close("Collection");
			xml.close("VTKFile");

			output.flush();

			if (!output) {
				throw std::runtime_error(
					"Failed to write PVD collection file: " +
					collection_path.string()
				);
			}
		}


	private:
		struct Frame {
			double time;
			std::filesystem::path file;
		};

		std::filesystem::path output_directory;
		std::filesystem::path collection_path;
		std::string base_name;
		std::vector<Frame> frames;
	};

} // namespace april