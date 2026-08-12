#pragma once

#include <concepts>
#include <type_traits>


namespace april {
    struct NoParticleAttributes {};
}

namespace april::particle {

    template <typename T>
    concept IsParticleAttributes =
        std::default_initializable<T> &&
        std::is_trivially_copyable_v<T> &&
        std::is_trivially_destructible_v<T> &&
        std::is_standard_layout_v<T> &&
        (!std::is_polymorphic_v<T>);


    // template struct used to tell the environment what user data will be used
    template<typename Data = NoParticleAttributes>
    struct ParticleAttributes { using particle_attributes_t = Data; };
    
}

namespace april {
    template<typename Data = NoParticleAttributes>
     inline constexpr particle::ParticleAttributes<Data> particle_attributes {};
}



