#pragma once
#include <filesystem>
#include <format>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <iostream>

#include "april/io/monitor.h"

namespace april::io {


	class TerminalOutput final : public Monitor {
	public:
		explicit TerminalOutput(const size_t write_frequency = 1): Monitor(write_frequency) {}

		void record(const size_t step, double, const Particles& particles) {
			std::cerr << "step: " << step <<  "\n";
			for (const auto & p : particles) {
				std::cerr << p.to_string() << "\n";
			}
		}
	};



	class BinaryOutput final : public Monitor {
	public:
		explicit BinaryOutput(const size_t write_frequency, std::string dir = "output", std::string base_name = "output"):
			Monitor(write_frequency), base_name(std::move(base_name)), dir(std::move(dir)) {}

		void record(size_t step, double, const Particles& particles) const {
			namespace fs = std::filesystem;

			fs::create_directories(dir);
			const std::string filename = std::format("{}_{:05}.bin", base_name, step);
			const fs::path full_path = fs::path(dir) / filename;

			std::ofstream out(full_path, std::ios::binary);
			if (!out) throw std::runtime_error("Failed to create output file: " + full_path.string());

			// Write header
			out.write(magic, sizeof(magic));							// 4 bytes
			write_binary(out, version);                             // 4 bytes
			write_binary(out, static_cast<uint64_t>(step));				// 8 bytes
			write_binary(out, static_cast<uint64_t>(particles.size()));	// 8 bytes
			write_binary(out, format_flags);                        // 4 bytes

			for (const auto& p : particles) {
				// Write position as 3 floats
				write_binary(out, static_cast<float>(p.position.x));
				write_binary(out, static_cast<float>(p.position.y));
				write_binary(out, static_cast<float>(p.position.z));

				write_binary(out, static_cast<uint32_t>(p.type));
				write_binary(out, static_cast<uint32_t>(p.id));
				write_binary(out, static_cast<uint8_t>(p.state));
			}
		}

		template<typename T> static void write_binary(std::ofstream& out, const T& value) {
			out.write(reinterpret_cast<const char*>(&value), sizeof(T));
		}

		std::string base_name;
		std::string dir;

		static constexpr char magic[4] = { 'P', 'A', 'R', 'T' };
		static constexpr uint32_t version = 1;
		static constexpr uint32_t format_flags = 0;
	};

} // namespace april::core