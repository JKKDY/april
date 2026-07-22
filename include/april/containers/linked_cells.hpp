#pragma once

#include "april/containers/layout/layout.hpp"

#include "april/containers/linked_cells/lc_aos.hpp"
#include "april/containers/linked_cells/lc_soa.hpp"
#include "april/containers/linked_cells/lc_aosoa.hpp"

namespace april {

    // user facing tag based instantiation
    template<typename Layout = Layout::AoSoA<>> class LinkedCells;


    template<>
    class LinkedCells<Layout::AoS> : public container::LinkedCellsAoS
    {};

    template<>
    class LinkedCells<Layout::SoA> : public container::LinkedCellsSoA
    {};

    template<uint8_t ChunkSize>
    class LinkedCells<Layout::AoSoA<ChunkSize>> : public container::LinkedCellsAoSoA<ChunkSize>
    {};

}




