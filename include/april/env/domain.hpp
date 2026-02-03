#pragma once

#include <optional>
#include "april/base/types.hpp"

namespace april::env {
	struct Domain {
		Domain() = default;
		Domain(const vec3d & originIn, const vec3d & extentIn): origin(originIn), extent(extentIn) {}

		std::optional<vec3d> origin;
		std::optional<vec3d> extent;

		[[nodiscard]] std::optional<double> volume() const noexcept {
			if (extent.has_value()) {
				return extent->x * extent->y * extent->z;
			}
			return std::nullopt;
		}

		[[nodiscard]] std::optional<vec3d> min_corner() const {
			if (origin && extent) {
				const vec3d far = origin.value() + extent.value();
				vec3d min;
				min.x = std::min(origin->x, far.x);
				min.y = std::min(origin->y, far.y);
				min.z = std::min(origin->z, far.z);
				return min;
			}

			return std::nullopt;
		}

		[[nodiscard]] std::optional<vec3d> max_corner() const {
			if (origin && extent) {
				const vec3d far = origin.value() + extent.value();
				vec3d max;
				max.x = std::max(origin->x, far.x);
				max.y = std::max(origin->y, far.y);
				max.z = std::max(origin->z, far.z);
				return max;
			}

			return std::nullopt;
		}

		[[nodiscard]] static Domain from_center_and_size(const vec3d& center, const vec3d& size) {
			const vec3d origin = center - (size / 2.0);
			return Domain{origin, size};
		}

		[[nodiscard]] static Domain from_min_max(const vec3d& min_corner, const vec3d& max_corner) {
			const vec3d origin = min_corner;
			const vec3d extent = max_corner - min_corner;
			return Domain{origin, extent};
		}
	};


	struct Box {
		Box() = default;

		[[nodiscard]] static Box from_domain(const Domain & domain) {
			if (domain.min_corner() && domain.max_corner()) {
				return {domain.min_corner().value(), domain.max_corner().value()};
			}
			throw std::logic_error("Domain not fully initialized. Got: "
				  "origin set: " + std::to_string(domain.origin.has_value()) +
				  "extent set: " + std::to_string(domain.extent.has_value())
			);
		}

		Box(const vec3d& min_corner, const vec3d& max_corner):
			min(min_corner), max(max_corner), extent(max-min) {

			AP_ASSERT(min_corner.x <= max_corner.x, "min_corner.x (" + std::to_string(min_corner.x) + ") "
				"is not <= than max_corner.x (" + std::to_string(max_corner.x)+")");
			AP_ASSERT(min_corner.y <= max_corner.y, "min_corner.y (" + std::to_string(min_corner.y) + ") "
				"is not <= than max_corner.y (" + std::to_string(max_corner.y)+")");
			AP_ASSERT(min_corner.z <= max_corner.z, "min_corner.z (" + std::to_string(min_corner.z) + ") "
				"is not <= than max_corner.z (" + std::to_string(max_corner.z)+")");
		}


		[[nodiscard]] bool contains(const vec3d & p) const noexcept{
			return (p.x >= min.x && p.x <= max.x) &&
				  (p.y >= min.y && p.y <= max.y) &&
				  (p.z >= min.z && p.z <= max.z);
		}

		[[nodiscard]] std::optional<Box> intersection(const Box& other) const {
			// The intersection of two boxes is the max of the mins and the min of the maxes.
			const double ix_min = std::max(min.x, other.min.x);
			const double iy_min = std::max(min.y, other.min.y);
			const double iz_min = std::max(min.z, other.min.z);

			const double ix_max = std::min(max.x, other.max.x);
			const double iy_max = std::min(max.y, other.max.y);
			const double iz_max = std::min(max.z, other.max.z);

			// If the new min is greater than the new max on ANY axis, there is no intersection.
			if (ix_min <= ix_max && iy_min <= iy_max && iz_min <= iz_max) {
				return Box(
					vec3d{ix_min, iy_min, iz_min},
					vec3d{ix_max, iy_max, iz_max}
					);
			}

			return std::nullopt;
		}

		[[nodiscard]] double volume() const noexcept {
			return extent.x * extent.y * extent.z;
		}

		const vec3d min;
		const vec3d max;
		const vec3d extent;
	};
}
