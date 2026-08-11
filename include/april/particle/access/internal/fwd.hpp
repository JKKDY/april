#pragma once
#include "april/particle/access/policy.hpp"
#include "april/particle/properties.hpp"
#include "april/particle/attributes.hpp"

namespace april::particle::internal {

    template <
        ParticleField ReadMask,
        ParticleField WriteMask,
        IsParticleAttributes Attributes,
        typename Source
    >
    struct PackedParticleRef;

    template <
       ParticleField ReadMask,
       ParticleField WriteMask,
       IsParticleAttributes Attributes,
       MaskPolicy MaskingPolicy
    >
    struct PackedParticleBuffer;

    template <
        ParticleField ReadMask,
        ParticleField WriteMask,
        IsParticleAttributes Attributes,
        MaskPolicy MaskingPolicy
    > struct PackedBufferView;
}