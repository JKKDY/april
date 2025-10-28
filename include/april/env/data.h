#pragma once


#include <vector>
#include <unordered_set>

#include "april/common.h"
#include "april/env/domain.h"
#include "april/env/particle.h"
#include "april/forces/force.h"
#include "april/boundaries/boundary.h"
#include "april/shared/pack_storage.h"

namespace april::env {
    inline const auto EXTENT_NOT_SET = vec3(std::numeric_limits<double>::max());
    inline const auto ORIGIN_NOT_SET = vec3(std::numeric_limits<double>::max());

    template<force::IsForcePack FPack,
            boundary::IsBoundaryPack BPack,
            controller::IsControllerPack CPack,
            field::IsFieldPack FFPack>
    class Environment;

}

namespace april::env::internal {

    // A data transfer object (DTO) to easily export environment data while preventing the user from directly
    // touching the environment internals

    template<class ForceVariant, class BoundaryVariant, class ControllerStorage, class FieldStorage>
    struct EnvironmentData final {
        Domain domain = {ORIGIN_NOT_SET, EXTENT_NOT_SET};

        // TODO document which takes precedent
        vec3 margin_abs = {0, 0, 0};
        vec3 margin_fac = {0.5, 0.5, 0.5}; // 50 % margin on each side by default

        std::unordered_set<env::ParticleID> user_particle_ids;
        std::unordered_set<env::ParticleType> user_particle_types;

        std::vector<env::Particle> particles;
        std::vector<force::internal::InteractionInfo<ForceVariant>> interactions {};
        std::array<BoundaryVariant, 6>  boundaries;

        ControllerStorage controllers;
        FieldStorage fields;
    };

    template<
            force::IsForcePack FPack,
            boundary::IsBoundaryPack BPack,
            controller::IsControllerPack CPack,
            field::IsFieldPack FFPack
        >
        auto get_env_data(const Environment<FPack, BPack, CPack, FFPack>& env) {
        return env.data;
    }


}