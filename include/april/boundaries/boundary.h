#pragma once

#include <cstdint>

#include "april/env/particle.h"

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

	inline int face_to_int(Face f) noexcept {
		return static_cast<int>(f);
	}

	inline uint8_t axis_of_face(const Face f) {
		return face_to_int(f) / 2;
	}

	inline bool face_sign_pos(const Face f) {
		return (face_to_int(f) & 1) != 0;
	}

	inline std::pair<uint8_t, uint8_t> non_face_axis(const Face f) {
		switch (axis_of_face(f)) {
		case 0: return {1,2};
		case 1: return {0,2};
		case 2: return {1,0};
		default: throw std::logic_error("got a different axis than {0,1,2}. This should never happen!?");
		}
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


