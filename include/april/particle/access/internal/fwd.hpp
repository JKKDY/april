#pragma once
#include "april/particle/access/policy.hpp"

namespace april::particle::internal {

    template <
        ParticleField ReadMask,
        ParticleField WriteMask,
        IsParticleAttributes Attributes
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