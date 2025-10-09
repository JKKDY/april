#pragma once

#include <cstdint>

#include "april/env/particle.h"
#include "april/env/domain.h"

namespace april::boundary {

	enum class Face : uint8_t {
		XMinus = 0, XPlus = 1,
		YMinus = 2, YPlus = 3,
		ZMinus = 4, ZPlus = 5,
	};

	constexpr std::array faces = {
		Face::XMinus, Face::XPlus,
		Face::YMinus, Face::YPlus,
		Face::ZMinus, Face::ZPlus
	};

	constexpr int to_int(Face f) noexcept {
		return static_cast<int>(f);
	}

	struct Topology {
		// Thickness of the boundary region adjacent to this face.
		//  > 0 : region lies inside the domain (e.g. reflective, repulsive)
		//  < 0 : region lies outside the domain (e.g. absorbing, teleporting)
		double boundary_thickness;

		// If true, this boundary couples its axis to the opposite face
		// (e.g. periodic boundaries: X- and X+ must both be periodic).
		bool couples_axis;

		// If true, this boundary changes iteration behaviour in the container
		// (e.g. periodic: requires min-image / ghost cells).
		// Otherwise, only particle dynamics are affected.
		bool force_wrap = false;
	};


	class Boundary {
	public:
		Boundary(const double thickness, const bool couples_axis, const bool force_wrap):
			topology(thickness, couples_axis, force_wrap)
		{}

		void dispatch_apply(this auto&& self, env::internal::Particle & particle) noexcept {
			static_assert(
			   requires { { self.apply(particle) } -> std::same_as<void>; },
			   "BoundaryCondition subclass must implement: void dispatch_apply(particle)"
			);

			self.apply(particle);
		}

		const Topology topology;
	};



	// define boundary concept
	template <class BC>
	concept IsBoundary = std::derived_from<BC, Boundary>;

	// define boundary pack
	template<IsBoundary... BCs>
	struct BoundaryPack {};

	template<class... BCs>
	requires (IsBoundary<BCs> && ...)
	inline constexpr BoundaryPack<BCs...> boundaries {};


	namespace internal {
		// define boundary variant concept
		template<typename T>
		struct is_boundary_variant : std::false_type {};

		// Accepts: std::variant<std::monostate, BCs...>
		template<typename... BCs>
		struct is_boundary_variant<std::variant<std::monostate, BCs...>>
		: std::bool_constant<(IsBoundary<BCs> && ...)> {};

		template<typename T>
		concept BoundaryVariant = is_boundary_variant<T>::value;
	}
}


