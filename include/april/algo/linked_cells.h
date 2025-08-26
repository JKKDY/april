#pragma once

#include "april/algo/algorithm.h"
#include "april/utils/set.hpp"
#include "april/env/particle.h"



namespace april::algo {
	namespace impl {
		class LinkedCells;
	}

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
		using impl = impl::LinkedCells;
		double cell_size_hint;
	};


	namespace impl {
		class LinkedCells final : public Algorithm<algo::LinkedCells> {
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
			explicit LinkedCells(const algo::LinkedCells & config);

			void build(const std::vector<Particle>& particles);
			void calculate_forces();

			[[nodiscard]] Particle & get_particle_by_id(ParticleID id);
			[[nodiscard]] ParticleID id_start() const;
			[[nodiscard]] ParticleID id_end() const;

			[[nodiscard]] Particle & get_particle_by_index(size_t index) noexcept;
			[[nodiscard]] size_t index_start() const;
			[[nodiscard]] size_t index_end() const;

			[[nodiscard]] size_t particle_count() const;
		private:
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
			std::vector<Particle> particles;
		};
	}
}