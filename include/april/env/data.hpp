#pragma once


#include <vector>
#include <unordered_set>
#include <array>

#include "april/base/types.hpp"
#include "april/env/domain.hpp"
#include "april/forces/force.hpp"
#include "april/boundaries/boundary.hpp"
#include "april/controllers/controller.hpp"
#include "april/fields/field.hpp"
#include "april/particle/generators.hpp"

namespace april::env {
    template<force::IsForcePack FPack,
            boundary::IsBoundaryPack BPack,
            controller::IsControllerPack CPack,
            field::IsFieldPack FFPack,
    	    IsUserData ParticleData>
    class Environment;

    struct ParticleCuboid;
    struct ParticleSphere;
}

namespace april::env::internal {

    // A data transfer object (DTO) to easily export environment data while preventing the user from directly
    // touching the environment internals

    // non templated base class
    struct EnvironmentCommonData  {
        Domain domain;

        // if both are specified the domain is chosen such that it satisfies both margin requirements
        vec3 margin_abs = {0, 0, 0};
        vec3 margin_fac = {0.5, 0.5, 0.5}; // 50 % margin on each side by default

        std::unordered_set<ParticleID> user_particle_ids;
        std::unordered_set<ParticleType> user_particle_types;

        std::vector<Particle> particles;
    };

    // templated extension
    template<class ForceVariant, class BoundaryVariant, class ControllerStorage, class FieldStorage>
    struct EnvironmentData final : EnvironmentCommonData{
        std::vector<force::internal::TypeInteraction<ForceVariant>> type_interactions {};
        std::vector<force::internal::IdInteraction<ForceVariant>> id_interactions {};
        std::array<BoundaryVariant, 6>  boundaries;

        ControllerStorage controllers;
        FieldStorage fields;
    };

    // friend function of environment to access the environment data
    template<
    force::IsForcePack FPack,
    boundary::IsBoundaryPack BPack,
    controller::IsControllerPack CPack,
    field::IsFieldPack FFPack,
    IsUserData ParticleData
    >
    auto get_env_data(const Environment<FPack, BPack, CPack, FFPack, ParticleData>& env) {
        return env.data;
    }
}