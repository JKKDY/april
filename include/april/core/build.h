#pragma once


#include <vector>
#include "april/core/system.h"
#include "april/env/particle.h"
#include "april/env/environment.h"
#include "april/containers/container.h"
#include "april/env/domain.h"


namespace april::core {
	namespace internal {
		struct InteractionParams {
			bool pair_contains_types;
			std::pair<int,int> key_pair;
		};

		// swap origin and extent components such that for all axis i: origin[i] < extent[i]
		env::Domain clean_domain(
			const env::Domain & domain
		);

		// calculate the minimal bounding box that contains all particles
		env::Box calculate_bounding_box(
			const std::vector<env::Particle>& particles
		);

		// check if user set domain parameters are ok (e.g. if it contains all particles)
		void validate_domain_params(
			const env::Domain& domain,
			const env::Box& bbox
		);

		// check if particle parameters are ok (e.g. no duplicate ids, every particle type has a force specified)
		void validate_particle_params(
			const std::vector<env::Particle> & particles,
			std::vector<InteractionParams> interactions,
			const std::unordered_set<env::ParticleID>& usr_particle_ids,
			const std::unordered_set<env::ParticleType>& usr_particle_types
		);

		// map user set particle ids & types to dense internal ids & types and return mappings
		UserToInternalMappings map_ids_and_types_to_internal(
			std::vector<env::Particle>& particles,
			const std::vector<InteractionParams>& interactions,
			std::unordered_set<env::ParticleID>& usr_particle_ids,
			std::unordered_set<env::ParticleType>& usr_particle_types
		);

		// given particle bounding box and user set parameters, calculate the final simulation domain
		env::Box finalize_environment_domain(
			const env::Box& bbox,
			const env::Domain& usr_domain,
			const vec3 & margin_abs,
			const vec3 & margin_fac
		);

		// build internal particle representation from user data
		std::vector<env::internal::Particle> build_particles(
			const std::vector<env::Particle>& particle_infos,
			const UserToInternalMappings& mapping
		);
	}

	template <container::IsContDecl C, class FPack, class BPack, class CPack, class FFPack>
	System<C, env::Environment<FPack, BPack, CPack, FFPack>> build_system(
		env::Environment<FPack, BPack, CPack, FFPack> & environment,
		const C& container,
		UserToInternalMappings* particle_mappings
	) {
		using EnvT = env::Environment<FPack, BPack, CPack, FFPack>;
		using BoundaryTable = boundary::internal::BoundaryTable<typename EnvT::boundary_variant_t>;
		using namespace internal;

		auto env = env::internal::get_env_data(environment);
		const env::Box bbox = calculate_bounding_box(env.particles);
		const env::Domain domain_cleaned = clean_domain(env.domain);

		std::vector<InteractionParams> interactions(env.interactions.size());
		for (size_t i = 0; i < interactions.size(); i++) {
			interactions[i].key_pair = env.interactions[i].key_pair;
			interactions[i].pair_contains_types = env.interactions[i].pair_contains_types;
		}

		validate_domain_params(
			domain_cleaned,
			bbox);

		validate_particle_params(
			env.particles,
			interactions,
			env.usr_particle_ids,
			env.usr_particle_types);

		const UserToInternalMappings mapping = map_ids_and_types_to_internal(
			env.particles,
			interactions,
			env.usr_particle_ids,
			env.usr_particle_types
		);

		const env::Box box = finalize_environment_domain(bbox, domain_cleaned, env.margin_abs, env.margin_fac);
		const env::Domain domain = {box.min, box.extent};
		const std::vector<env::internal::Particle> particles = build_particles(env.particles, mapping);

		if (particle_mappings) {
			particle_mappings->usr_ids_to_impl_ids = mapping.usr_ids_to_impl_ids;
			particle_mappings->usr_types_to_impl_types = mapping.usr_types_to_impl_types;
		}

		BoundaryTable boundaries (env.boundaries, domain);

		std::vector<boundary::Topology> topologies;
		for (boundary::Face face : boundary::all_faces) {
			topologies.push_back(boundaries.get_boundary(face).topology);
		}

		// TODO validate topologies

		container::internal::ContainerFlags container_flags = {};
		for (const boundary::Face face : boundary::all_faces) {
			if (topologies[face_to_int(face)].force_wrap) {
				switch (axis_of_face(face)) {
				case 0: container_flags.periodic_x = true; break;
				case 1: container_flags.periodic_y = true; break;
				case 2: container_flags.periodic_z = true; break;
				default: std::unreachable();
				}
			}
		}

		return System<C, env::Environment<FPack, BPack, CPack, FFPack>> (
			container,
			container_flags,
			domain,
			particles,
			boundaries,
			env.controllers,
			env.fields,
			mapping.usr_types_to_impl_types,
			mapping.usr_ids_to_impl_ids,
			env.interactions
		);
	}
}