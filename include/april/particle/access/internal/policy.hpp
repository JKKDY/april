#pragma once

namespace april::particle::internal {
    enum class MaskPolicy {
        Disabled,
        Enabled
    };


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