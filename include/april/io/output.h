#pragma once
#include <concepts>
#include <filesystem>
#include <format>
#include <fstream>
#include <string>
#include <vector>

#include "april/env/particle.h"

namespace april::io {

	template <typename T> concept IsOutputWriter = requires(T t, double time,
		const std::vector<env::impl::Particle>& particles) {
		{ t.write(time, particles) } -> std::same_as<void>;
	};

	class OutputWriter {
	public:
		explicit OutputWriter(const size_t write_frequency)
			: write_frequency(write_frequency) {}

		void write(this auto&& self, size_t step, const std::vector<env::impl::Particle>& particles) {
			self.write_output(step, particles); // works if 'self' has a write(t, particles) method
		}

		const size_t write_frequency;
	};

	class NullOutput final : public OutputWriter {
		void write_output(size_t, const std::vector<env::impl::Particle>&) {}
	};

	class BinaryOutput final : public OutputWriter {
	public:
		explicit BinaryOutput(const size_t write_frequency, const std::string& dir = "output",const std::string& base_name = "output"):
			OutputWriter(write_frequency), base_name(base_name), dir(dir) {}

	private:
		void write_output(size_t step, const std::vector<env::impl::Particle>& particles) const {
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

		const std::string base_name;
		const std::string dir;

		static constexpr char magic[4] = { 'P', 'A', 'R', 'T' };
		static constexpr uint32_t version = 1;
		static constexpr uint32_t format_flags = 0;
	};

} // namespace april::core