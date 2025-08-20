#pragma once

#include "april/algo/container.h"
#include "april/utils/set.hpp"
#include "april/env/particle.h"

namespace april::core::impl {
	class LinkedCells;
}

namespace april::core {


	template<typename T> concept IsParticleSet = requires(T set, size_t id) {
		// Construction with maxId
		{ T(std::declval<size_t>()) };

		// Insert method
		{ set.insert(id) } -> std::same_as<void>;

		// Erase method
		{ set.erase(id) } -> std::same_as<void>;

		// Contains method
		{ set.contains(id) } -> std::same_as<bool>;

		// Iteration methods
		{ set.begin() } -> std::input_iterator;
		{ set.end() } -> std::input_iterator;

		// Size method
		{ set.size() } -> std::same_as<size_t>;
	};


	struct LinkedCells {
		using Container = impl::LinkedCells;
		double cell_size_hint;
	};


	namespace impl {
		class LinkedCells final : public Container {
			using Particle = env::impl::Particle;

			struct Cell {
				using ParticleSet = utils::IndexSet<env::impl::ParticleID>;
				static_assert(IsParticleSet<ParticleSet>, "ParticleSet must implement IsParticleSet interface");

				ParticleSet particles;
				uint3 idx;
				unsigned int id;
			};

			struct CellPair {
				Cell& first;
				Cell& second;
			};

		public:
			explicit LinkedCells(double cell_size = -1);

			void build() override;
			void calculate_forces() override;

		private:
			double grid_constant;
			void build_cells();
			void build_cell_pairs();

			Cell & get_cell(const vec3& position) noexcept;
			Cell & get_cell(const uint3 & idx);

			vec3 cell_size;
			vec3 inv_cell_size;
			uint3 cell_count{};
			Cell outside_cell;
			std::vector<Cell> cells;
			std::vector<CellPair> cell_pairs;
		};
	}
}