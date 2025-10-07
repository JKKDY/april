#pragma once

#include "april/common.h"

namespace april::env {
	struct Domain {
		vec3 extent = {};
		vec3 origin = {};

		[[nodiscard]] bool contains(const vec3 & p) const noexcept{
			const vec3 max = extent + origin;
			return (p.x >= origin.x && p.x <= max.x) &&
				  (p.y >= origin.y && p.y <= max.y) &&
				  (p.z >= origin.z && p.z <= max.z);
		}

		[[nodiscard]] double volume() const noexcept {
			return extent.x * extent.y * extent.z;
		}
	};


}
// move boundary to domain/