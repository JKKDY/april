#pragma once

#include "april/common.h"

namespace april::env {
	struct Domain {
		vec3 extent = {};
		vec3 origin = {};

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
		min(domain.min_corner()), max(domain.max_corner()), extent(max-min){
			AP_ASSERT(vec3::all(domain.origin,
				[=](auto x) {return x <= std::numeric_limits<double>::max() / 2 &&  x >= std::numeric_limits<double>::lowest() / 2;}),
				"all modulus components in origin must be within range [min/2, max/2] to avoid overflow");
			AP_ASSERT(vec3::all(domain.extent,
				[=](auto x) {return x <= std::numeric_limits<double>::max() / 2 &&  x >= std::numeric_limits<double>::lowest() / 2;}),
				"all modulus components in origin must be within range [min/2, max/2] to avoid overflow");
		}

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
