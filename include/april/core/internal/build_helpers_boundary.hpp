#pragma once


#include "april/boundaries/boundary.hpp"
#include "april/containers/container.hpp"


namespace april::core::internal {

	// Extracts the topology data for each of the 6 domain faces
	template<class BoundaryTable>
		auto extract_topologies(const BoundaryTable  & boundaries) {
		std::vector<boundary::Topology> topologies;
		for (DomainFace face : all_faces) {
			topologies.push_back(boundaries[face].topology);
		}
		return topologies;
	}

	// Populates unconfigured faces with a default OpenBoundary
	template<boundary::internal::IsBoundaryVariant BV>
		auto set_default_boundaries(std::array<BV, 6>  & boundaries) {
		for (auto & v : boundaries)
			if (std::holds_alternative<boundary::internal::BoundarySentinel>(v))
				v.template emplace<OpenBoundary>(); // default-construct Open boundary
	}


	// Ensures coupled boundaries (e.g. Periodic) are applied symmetrically to an axis
	inline void validate_topologies(const std::vector<boundary::Topology> & topologies) {
		for (int axis = 0; axis < 3; ++axis) {
			const auto & face_min = topologies[axis * 2];
			const auto & face_plus = topologies[axis * 2 + 1];

			// If either face requires coupling (Periodic, Teleport, etc.), they must both agree
			if (face_min.couples_axis || face_plus.couples_axis) {
				if (face_min.couples_axis != face_plus.couples_axis ||
					face_min.force_wrap != face_plus.force_wrap) {

					const std::string axis_name = (axis == 0) ? "X" : (axis == 1) ? "Y" : "Z";
					throw std::invalid_argument(std::format(
					   "Asymmetric boundary configuration on {} axis. Coupled boundaries "
					   "(like Periodic) must be applied to both - and + faces.",
					   axis_name
					));
				}
			}
		}
	}


	// Maps boundary topology requirements to backend container configuration flags
	inline container::ContainerFlags set_container_flags(const std::vector<boundary::Topology>& topologies) {
		container::ContainerFlags container_flags = {};
		for (const DomainFace face : all_faces) {
			if (topologies[boundary::face_to_int(face)].force_wrap) {
				switch (boundary::axis_of_face(face)) {
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















