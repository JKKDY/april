#pragma once
#include <filesystem>
#include <format>
#include <fstream>
#include <string>
#include <utility>


#include "april/monitors/monitor.h"


namespace april::monitor {

	class BinaryOutput final : public Monitor {
	public:
		static constexpr env::FieldMask fields = to_field_mask(env::Field::all);
		explicit BinaryOutput(
			const shared::Trigger & trigger,
			std::string dir = "output",
			std::string base_name = "output")
		:
			Monitor(trigger), base_name(std::move(base_name)), dir(std::move(dir)) {}

		template<class S>
		void record(const core::SystemContext<S> & sys) const {
			namespace fs = std::filesystem;

			size_t start_idx = sys.index_start();
			size_t end_idx = sys.index_end();

			fs::create_directories(dir);
			const std::string filename = std::format("{}_{:05}.bin", base_name, sys.step());
			const fs::path full_path = fs::path(dir) / filename;

			std::ofstream out(full_path, std::ios::binary);
			if (!out) throw std::runtime_error("Failed to create output file: " + full_path.string());

			// Write header
			out.write(magic, sizeof(magic));				// 4 bytes
			write_binary(out, version);                 // 4 bytes
			write_binary(out, sys.step());				// 8 bytes
			write_binary(out, end_idx - start_idx);		// 8 bytes
			write_binary(out, format_flags);            // 4 bytes

			for (size_t i = start_idx; i < end_idx; ++i) {
				auto p = sys.template get_particle_by_index<fields>(i);
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

	private:
		std::string base_name;
		std::string dir;

		static constexpr char magic[4] = { 'P', 'A', 'R', 'T' };
		static constexpr uint32_t version = 1;
		static constexpr uint32_t format_flags = 0;
	};

} // namespace april::core