#pragma once

#include "april/containers/layout.hpp"

#include "april/containers/direct_sum/ds_aos.hpp"
#include "april/containers/direct_sum/ds_soa.hpp"
#include "april/containers/direct_sum/ds_aosoa.hpp"

namespace april {

    // User-facing tag-based instantiation with AoSoA as the default
    template<typename Layout = Layout::AoSoA<>>
    class DirectSum;

    template<>
    class DirectSum<Layout::AoS> : public container::DirectSumAoS
    {};

    template<>
    class DirectSum<Layout::SoA> : public container::DirectSumSoA
    {};

    template<>
    class DirectSum<Layout::AoSoA<>> : public container::DirectSumAoSoA
    {};

}