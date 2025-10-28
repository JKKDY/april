#pragma once

#include <optional>
#include "april/common.h"

namespace april::env {
	inline const auto EXTENT_NOT_SET = vec3(std::numeric_limits<double>::max());
	inline const auto ORIGIN_NOT_SET = vec3(std::numeric_limits<double>::max());

	struct Domain {
		Domain() = default;
		Domain(const vec3 & originIn, const vec3 & extentIn): origin(originIn), extent(extentIn) {}

		std::optional<vec3> origin;
		std::optional<vec3> extent;

		[[nodiscard]] std::optional<double> volume() const noexcept {
			if (extent.has_value()) {
				return extent->x * extent->y * extent->z;
			}
			return std::nullopt;
		}

		[[nodiscard]] std::optional<vec3> min_corner() const {
			if (origin && extent) {
				const vec3 far = origin.value() + extent.value();
				vec3 min;
				min.x = std::min(origin->x, far.x);
				min.y = std::min(origin->y, far.y);
				min.z = std::min(origin->z, far.z);
				return min;
			}

			return std::nullopt;
		}

		[[nodiscard]] std::optional<vec3> max_corner() const {
			if (origin && extent) {
				const vec3 far = origin.value() + extent.value();
				vec3 max;
				max.x = std::max(origin->x, far.x);
				max.y = std::max(origin->y, far.y);
				max.z = std::max(origin->z, far.z);
				return max;
			}

			return std::nullopt;
		}

		[[nodiscard]] static Domain from_center_and_size(const vec3& center, const vec3& size) {
			const vec3 origin = center - (size / 2.0);
			return Domain{origin, size};
		}

		[[nodiscard]] static Domain from_min_max(const vec3& min_corner, const vec3& max_corner) {
			const vec3 origin = min_corner;
			const vec3 extent = max_corner - min_corner;
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

			// 	AP_ASSERT(vec3::all(domain.origin,
			// 		[=](auto x) {return x <= std::numeric_limits<double>::max() / 2 &&  x >= std::numeric_limits<double>::lowest() / 2;}),
			// 		"all modulus components in origin must be within range [min/2, max/2] to avoid overflow");
			// 	AP_ASSERT(vec3::all(domain.extent,
			// 		[=](auto x) {return x <= std::numeric_limits<double>::max() / 2 &&  x >= std::numeric_limits<double>::lowest() / 2;}),
			// 		"all modulus components in origin must be within range [min/2, max/2] to avoid overflow");
			// }


			Box(const vec3& min_corner, const vec3& max_corner):
				min(min_corner), max(max_corner), extent(max-min) {

				AP_ASSERT(min_corner.x <= max_corner.x, "min_corner.x (" + std::to_string(min_corner.x) + ") "
					"is not <= than max_corner.x (" + std::to_string(max_corner.x)+")");
				AP_ASSERT(min_corner.y <= max_corner.y, "min_corner.y (" + std::to_string(min_corner.y) + ") "
					"is not <= than max_corner.y (" + std::to_string(max_corner.y)+")");
				AP_ASSERT(min_corner.z <= max_corner.z, "min_corner.z (" + std::to_string(min_corner.z) + ") "
					"is not <= than max_corner.z (" + std::to_string(max_corner.z)+")");
			}


			[[nodiscard]] bool contains(const vec3 & p) const noexcept{
				return (p.x >= min.x && p.x <= max.x) &&
					  (p.y >= min.y && p.y <= max.y) &&
					  (p.z >= min.z && p.z <= max.z);
			}

			const vec3 min;
			const vec3 max;
			const vec3 extent;
		};
}
