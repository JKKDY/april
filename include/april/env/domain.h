#pragma once

#include "april/common.h"

namespace april::env {
	struct Domain {
		Domain() = default;
		Domain(const vec3 & origin_, const vec3 & extent_): origin(origin_), extent(extent_) {}

		vec3 origin = {};
		vec3 extent = {};

		[[nodiscard]] double volume() const noexcept {
			return extent.x * extent.y * extent.z;
		}

		[[nodiscard]] vec3 min_corner() const {
			const vec3 far = origin + extent;

			vec3 min;
			min.x = std::min(origin.x, far.x);
			min.y = std::min(origin.y, far.y);
			min.z = std::min(origin.z, far.z);

			return min;
		}

		[[nodiscard]] vec3 max_corner() const {
			const vec3 far = origin + extent;

			vec3 max;
			max.x = std::max(origin.x, far.x);
			max.y = std::max(origin.y, far.y);
			max.z = std::max(origin.z, far.z);

			return max;
		}
	};

	struct Box {
		explicit Box(const Domain & domain):
		Box(domain.min_corner(), domain.max_corner()){
			AP_ASSERT(vec3::all(domain.origin,
				[=](auto x) {return x <= std::numeric_limits<double>::max() / 2 &&  x >= std::numeric_limits<double>::lowest() / 2;}),
				"all modulus components in origin must be within range [min/2, max/2] to avoid overflow");
			AP_ASSERT(vec3::all(domain.extent,
				[=](auto x) {return x <= std::numeric_limits<double>::max() / 2 &&  x >= std::numeric_limits<double>::lowest() / 2;}),
				"all modulus components in origin must be within range [min/2, max/2] to avoid overflow");
		}

		Box(const vec3& min_corner, const vec3& max_corner):
		min(min_corner), max(max_corner), extent(max-min){}


		[[nodiscard]] bool contains(const vec3 & p) const noexcept{
			return (p.x >= min.x && p.x <= max.x) &&
				  (p.y >= min.y && p.y <= max.y) &&
				  (p.z >= min.z && p.z <= max.z);
		}

		vec3 min;
		vec3 max;
		vec3 extent;
	};
}
