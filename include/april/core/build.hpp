#pragma once

#include <unordered_map>
#include <utility>

#include "april/core/system.hpp"
#include "april/core/environment.hpp"
#include "april/core/domain.hpp"
#include "april/core/internal/build_helpers_particle.hpp"
#include "april/core/internal/build_helpers_domain.hpp"
#include "april/core/internal/build_helpers_boundary.hpp"
#include "april/containers/container.hpp"


namespace april {
    // @brief build-time metadata for system
    struct BuildInfo {
        std::unordered_map<ParticleType, ParticleType> type_map; // user type -> dense system type
        std::unordered_map<ParticleID, ParticleID> id_map; // user id -> dense system id
        Domain particle_box; // minimal bounding box containing all particles
        Domain simulation_domain; // final domain including margins
    };


    // @brief constructs a system from an environment, container and execution config. Returns optional build information via build_info pointer
    template <class ContainerCfg, core::IsEnvironment Env, exec::IsExecutionConfig ExecCfg>
        requires container::IsContainerDecl<ContainerCfg, typename Env::traits, ExecCfg>
    auto build_system(
        const Env& environment,
        const ContainerCfg& container_config,
        const ExecCfg& execution_config,
        BuildInfo* build_info
    ) {
        using namespace april::core::internal;
        using BoundaryTable = Env::traits::boundary_table_t;
        using ForceTable = Env::traits::force_table_t;
        using ParticleAttributes = Env::traits::particle_attributes_t;

        // explicit type for IDE code completion
        using EnvData = EnvironmentData<
            typename Env::traits::force_variant_t,
            typename Env::traits::boundary_variant_t,
            typename Env::traits::controller_storage_t,
            typename Env::traits::field_storage_t>;

        // copy environment data
        EnvData env = get_env_data(environment);

        // validate & set simulation domain
        const core::Box particle_bbox = particle_bounding_box(env.particles);
        const core::Box simulation_box = determine_simulation_box(
            env.domain, particle_bbox, env.margin_abs, env.margin_fac);

        // get all interacting type and id pairs
        auto [type_pairs, id_pairs] =
            extract_interaction_parameters(env.type_interactions, env.id_interactions);

        // ensure unique ids and create user -> system mappings for ids and types
        assign_missing_particle_ids(env.particles, env.user_particle_ids);
        auto [type_map, id_map] = create_particle_mappings(
            env.particles,
            env.user_particle_types,
            env.user_particle_ids,
            type_pairs,
            id_pairs
        );

        // create particles
        auto particles = build_particles<ParticleAttributes>(env.particles, type_map, id_map);

        // create force table
        ForceTable forces(env.type_interactions, env.id_interactions, type_map, id_map);

        // if no boundary specified use a default (OpenBoundary)
        set_default_boundaries(env.boundaries);
        BoundaryTable boundaries(env.boundaries, simulation_box);
        auto topologies = extract_topologies(boundaries);
        validate_topologies(topologies);

        // fill build info if requested
        if (build_info) {
            build_info->type_map = type_map;
            build_info->id_map = id_map;
            build_info->particle_box = Domain(particle_bbox.min, particle_bbox.extent);
            build_info->simulation_domain = Domain(simulation_box.min, simulation_box.extent);
        }

        // build container
        using ContainerConfig = container::ContainerBuildConfig<ContainerCfg, ExecCfg, ParticleAttributes>;
        ContainerConfig container_build_config{
            .exec = execution_config,
            .config = container_config,
            .flags = set_container_flags(topologies),
            .hints = container::ContainerHints(),
            .interaction_map = forces.generate_interaction_map(),
            .domain = simulation_box
        };

        auto container = typename ContainerConfig::Container(container_build_config);

        // build system
        SystemBuildConfig<typename ContainerConfig::Container, typename Env::traits, ExecCfg> system_config{
            .container = std::move(container),
            .execution_config = execution_config,
            .particles = std::move(particles),
            .boundaries = std::move(boundaries),
            .interactions = std::move(forces),
            .controllers = std::move(env.controllers),
            .fields = std::move(env.fields)
        };

        return System(std::move(system_config));
    }


    // convenience overloads
    template <class Container, core::IsEnvironment EnvT>
    auto build_system(
        const EnvT& environment,
        const Container& container_config
    ) {
        return build_system(environment, container_config, ExecutionConfig(), nullptr);
    }

    template <class Container, core::IsEnvironment EnvT>
    auto build_system(
        const EnvT& environment,
        const Container& container_config,
        BuildInfo* build_info
    ) {
        return build_system(environment, container_config, ExecutionConfig(), build_info);
    }

    template <class Container, core::IsEnvironment EnvT>
    auto build_system(
        const EnvT& environment,
        const Container& container_config,
        const exec::IsExecutionConfig auto& execution_config
    ) {
        return build_system(environment, container_config, execution_config, nullptr);
    }
}
