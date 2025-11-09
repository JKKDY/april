#pragma once


#include "april/boundaries/boundary.h"
#include "april/containers/container.h"


namespace april::core::internal {
	// --- set boundaries
	template<class BoundaryTable>
		auto extract_topologies(const BoundaryTable  & boundaries) {
		std::vector<boundary::Topology> topologies;
		for (boundary::Face face : boundary::all_faces) {
			topologies.push_back(boundaries.get_boundary(face).topology);
		}
		return topologies;
	}


	template<boundary::internal::IsBoundaryVariant BV>
		auto set_default_boundaries(std::array<BV, 6>  & boundaries) {
		for (auto & v : boundaries)
			if (std::holds_alternative<boundary::internal::BoundarySentinel>(v))
				v.template emplace<boundary::Open>(); // default-construct Open boundary
	}


	inline void validate_topologies(const std::vector<boundary::Topology> &) {
		// TODO implement validate_topologies
	}


	inline container::internal::ContainerFlags set_container_flags(const std::vector<boundary::Topology>& topologies) {
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
		return container_flags;
	}
}
