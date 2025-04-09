#pragma once


#pragma once
#include <stdint.h>

#include "ankerl/unordered_dense.h"
#include "common.h"
#include "particle.h"

namespace april::env {

    namespace impl {
        struct PointerHash {
            size_t operator()(const void* ptr) const noexcept {
                return reinterpret_cast<uintptr_t>(ptr) >> 3; // Right shift to improve distribution
            }
        };



        struct GridCell {
            using Cell_T = ankerl::unordered_dense::set<Particle*, PointerHash>;

            enum Type {
                INNER = 0x001,              ///< cell inside the grid.
                BOUNDARY_RIGHT = 0x010,     ///< boundary cell on the right side.
                BOUNDARY_LEFT = 0x020,      ///< boundary cell on the left side.
                BOUNDARY_TOP = 0x040,       ///< boundary cell on the top side.
                BOUNDARY_BOTTOM = 0x080,    ///< boundary cell on the bottom side.
                BOUNDARY_FRONT = 0x100,     ///< boundary cell on the front side.
                BOUNDARY_BACK = 0x200,      ///< boundary cell on the back side.
                BOUNDARY = 0xFF0,           ///< boundary cell.
                OUTSIDE = 0x002,            ///< cell outside the grid.
                INSIDE = INNER | BOUNDARY,  ///< Cells considered inside the boundary
            };

            GridCell(const vec3& coord, const vec3& size, Type type, const int3& idx);


            bool operator==(const GridCell& other) const;

            std::string to_string() const;

            const Type type;   ///< The type of the grid cell.
            const vec3 origin; ///< The coordinates of the origin of the grid cell.
            const vec3 size;   ///< The size of the grid cell.
            const int3 idx;    ///< The index of the grid cell
            Cell_T particles;  ///< The particles in the grid cell.
        };


        struct CellPair {
            enum Periodicity {
                PERIODIC_NONE = 0x0,  ///< No periodicity.
                PERIODIC_X = 0x1,  ///< Periodicity in the x-direction.
                PERIODIC_Y = 0x2,  ///< Periodicity in the y-direction.
                PERIODIC_Z = 0x4,  ///< Periodicity in the z-direction.
            };

            CellPair(GridCell& cell1, GridCell& cell2, Periodicity periodicity);

            bool empty() const;
            std::string to_string() const;
            std::pair<int, int> id() const; //Retrieves the ids of the cells

            const Periodicity periodicity;
            GridCell& cell1;  ///< The first grid cell in the pair.
            GridCell& cell2;  ///< The second grid cell in the pair.
        };



        class ParticleGrid {
        public:
            ParticleGrid();
            void build(const Boundary& boundary, double grid_const, std::vector<Particle>& particles, bool build_blocks = false);


            GridCell& get_cell(const int3& idx);
            GridCell& get_cell(const Particle& particle);

            [[nodiscard]] const GridCell& get_cell(const int3& idx) const;
            [[nodiscard]] const GridCell& get_cell(const Particle& particle) const;
            [[nodiscard]] int3 what_cell(const vec3& position) const;
            [[nodiscard]] std::vector<int3> get_cell_indices() const;


            const std::vector<CellPair>& linked_cells();
            const std::vector<GridCell*>& boundary_cells();


            size_t particle_count() const;
            void update_cells(Particle* particle, const int3& old_cell, const int3& new_cell);

        private:

            void build_cells(const vec3& extent, double grid_constant, std::vector<Particle>& particles);
            void build_cell_pairs(const std::array<BoundaryRule, 6>& rules);

            ankerl::unordered_dense::map<int3, GridCell, Int3Hasher> cells{};  ///< A hash map storing the cells in the grid.
            std::vector<CellPair> cell_pairs{};                  ///< A vector of linked cell pairs.
            std::vector<GridCell*> border_cells;                 ///< A vector of cells at the domain boundary

            uint3 cell_count{};         ///< The number of cells in the grid along each dimension.
            vec3 cell_size{};           ///< The size of each grid cell.
            vec3 boundary_origin = {};  ///< The origin of the boundary.
        };




        //Enables bitwise OR operation for GridCell::Type.
        constexpr GridCell::Type operator|(const GridCell::Type lhs, const GridCell::Type rhs) {
            using T = std::underlying_type_t<GridCell::Type>;
            return static_cast<GridCell::Type>(static_cast<T>(lhs) | static_cast<T>(rhs));
        }

        //Enables bitwise OR operation for GridCell::Type.
        constexpr bool operator&(const GridCell::Type lhs, const GridCell::Type rhs) {
            using T = std::underlying_type_t<GridCell::Type>;
            return static_cast<T>(lhs) & static_cast<T>(rhs);
        }
        //Enables bitwise OR assignment for GridCell::Type.
        constexpr GridCell::Type& operator|=(GridCell::Type& lhs, const GridCell::Type rhs) {
            using T = std::underlying_type_t<GridCell::Type>;
            lhs = static_cast<GridCell::Type>(static_cast<T>(lhs) | static_cast<T>(rhs));
            return lhs;
        }

        // Enables bitwise OR operation for combining periodicity flags.
        inline CellPair::Periodicity operator|(const CellPair::Periodicity lhs, const CellPair::Periodicity rhs) {
            using T = std::underlying_type_t<CellPair::Periodicity>;
            return static_cast<CellPair::Periodicity>(static_cast<T>(lhs) | static_cast<T>(rhs));
        }
        // Enables bitwise OR assignment for periodicity flags.
        inline CellPair::Periodicity& operator|=(CellPair::Periodicity& lhs, CellPair::Periodicity rhs) {
            using T = std::underlying_type_t<CellPair::Periodicity>;
            lhs = static_cast<CellPair::Periodicity>(static_cast<T>(lhs) | static_cast<T>(rhs));
            return lhs;
        }
        // Enables bitwise AND operator for periodicity flags.
        inline CellPair::Periodicity operator&(const CellPair::Periodicity lhs, const CellPair::Periodicity rhs) {
            using T = std::underlying_type_t<CellPair::Periodicity>;
            return static_cast<CellPair::Periodicity>(static_cast<T>(lhs) & static_cast<T>(rhs));
        }
    } 
}  // namespace april::env