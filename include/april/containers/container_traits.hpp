#pragma once

#include "april/core/domain.hpp"
#include "april/particle/attributes.hpp"
#include "april/particle/record.hpp"
#include "april/interactions/interaction_table.hpp"
#include <vector>

#include "april/exec/config.hpp"

namespace april::container {
    struct ContainerFlags {
        bool periodic_x;			// domain is periodic along x-axis
        bool periodic_y;			// domain is periodic along y-axis
        bool periodic_z;			// domain is periodic along z-axis
        bool infinite_domain;		// particles outside of domain still interact normally (time complexity may go to O(n^2)
        bool particle_addable;		// particles can be added during run time
        bool particle_deletable;	// particles can be deleted during run time
    };

    struct ContainerHints {
        // TODO add regions that will be queried in the future so the container can keep track of particles better
        std::vector<ParticleID> interacting_particles;
        std::vector<core::Box> query_regions;
    };

    template<class ContainerCfg, class ExecutionCfg, particle::IsParticleAttributes Attributes>
    struct ContainerBuildConfig {
        using ContainerConfig = ContainerCfg;
        using ExecutionConfig = ExecutionCfg;
        using ParticleAttributes = Attributes;
        using ThreadExecutor = ExecutionConfig::ThreadExecutor;
        using Container = ContainerCfg::template impl<ContainerBuildConfig>;

        ExecutionConfig exec;
        ContainerConfig config;
        ContainerFlags flags {};
        ContainerHints hints {};
        interactions::internal::InteractionMap interaction_map {};
        core::Box domain {};
    };


    namespace internal {
		// check if T is a ContainerBuildConfig
		template<typename T>
		struct is_container_build_context : std::false_type {};

		template<typename Cfg, typename Exec, typename Attributes>
		struct is_container_build_context<
			ContainerBuildConfig<Cfg, Exec, Attributes>
		> : std::true_type {};
	}

	template<typename T>
	concept IsContainerBuildConfig =
		internal::is_container_build_context<std::remove_cvref_t<T>>::value;


	// forward declaration
	template<IsContainerBuildConfig BuildConfiguration>
	class Container;



	template<typename C>
	concept HasContainerOps = requires (
	    C c,
	    const C cc,
	    ParticleID id,
	    size_t index,
	    const core::Box& region,
	    const particle::ParticleRecord<typename C::ParticleAttributes>& p,
	    const std::vector<particle::ParticleRecord<typename C::ParticleAttributes>>& particles
	) {
		// minimal implemented interface (except get_field_ptr since that is not part of the public API)
	    { c.build(particles) };
	    { c.rebuild_structure() };

		{ cc.capacity() } -> std::convertible_to<size_t>;
	    { cc.particle_count() } -> std::convertible_to<size_t>;
		{ cc.min_id() } -> std::convertible_to<ParticleID>;
	    { cc.max_id() } -> std::convertible_to<ParticleID>;

	    { cc.id_to_index(id) } -> std::convertible_to<size_t>;
		{ cc.contains_id(id) } -> std::convertible_to<bool>;
		{ cc.index_is_valid(id) } -> std::convertible_to<bool>;

	    { c.collect_indices_in_region(region) } -> std::convertible_to<std::vector<size_t>>;

	    { c.template for_each_interaction_batch<ParallelPolicy::Serial>([](auto&&){}) };
		{ c.template for_each_topology_batch<ParallelPolicy::Serial>([](auto&&){}) };
	};


	template<typename C>
	concept IsContainer =
		requires {
			typename C::Config;
			typename C::ExecutionConfig;
			typename C::ParticleAttributes;
		}
		&& HasContainerOps<C> // has implemented minimal interface
		&& std::derived_from<C, // is derived from container base class
			Container<ContainerBuildConfig<typename C::Config, typename C::ExecutionConfig, typename C::ParticleAttributes>>
		>;


	template<typename ContainerDecl, typename Traits, typename ExecCfg>
	concept IsContainerDecl =
		core::internal::IsEnvironmentTraits<Traits> &&
		exec::IsExecutionConfig<ExecCfg> &&
		requires { // check if ContainerDecl has an impl member that is template-able on a build configuration
			typename ContainerBuildConfig<ContainerDecl, ExecCfg, typename Traits::particle_attributes_t>;

			typename ContainerDecl::template impl<
				ContainerBuildConfig<ContainerDecl, ExecCfg, typename Traits::particle_attributes_t>
			>;
		} &&
		IsContainer< // and has represents a valid container
			typename ContainerDecl::template impl<
				ContainerBuildConfig<ContainerDecl, ExecCfg, typename Traits::particle_attributes_t>
			>
		>;
}
