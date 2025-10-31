#pragma once


#include <vector>
#include <unordered_set>

#include "april/common.h"
#include "april/env/domain.h"
#include "april/env/particle.h"
#include "april/forces/force.h"
#include "april/boundaries/boundary.h"
#include "april/controllers/controller.h"
#include "april/fields/field.h"

namespace april::env {
    template<force::IsForcePack FPack,
            boundary::IsBoundaryPack BPack,
            controller::IsControllerPack CPack,
            field::IsFieldPack FFPack,
    	    IsParticleData ParticleData>
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
    IsParticleData ParticleData
    >
    auto get_env_data(const Environment<FPack, BPack, CPack, FFPack, ParticleData>& env) {
        return env.data;
    }

    void add_particle_impl(EnvironmentCommonData& data, const Particle& particle);
    std::vector<ParticleID> add_cuboid_particles_impl(EnvironmentCommonData& data, const ParticleCuboid& cuboid);
    std::vector<ParticleID> add_sphere_particles_impl(EnvironmentCommonData& data, const ParticleSphere& sphere);

}