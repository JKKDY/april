#pragma once

#include <cstdint>
#include <concepts>

#include "april/base/traits.hpp"
#include "april/particle/defs.hpp"
#include "april/particle/access.hpp"
#include "april/env/domain.hpp"

namespace april::boundary {

	struct Open;

	enum class Face : uint8_t {
		XMinus = 0, XPlus = 1,
		YMinus = 2, YPlus = 3,
		ZMinus = 4, ZPlus = 5,
	};

	const std::vector all_faces = {
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
		explicit Topology(
			const double thickness,
			const bool coupled_axis,
			const bool has_force_wrap,
			const bool change_particle_pos)
		:
			boundary_thickness(thickness),
			couples_axis(coupled_axis),
			force_wrap(has_force_wrap),
			may_change_particle_position(change_particle_pos) {}

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
		bool force_wrap;

		// If true, the boundary the containers register_particle_movement routine
		// will be called after each invocation of the boundary condition
		bool may_change_particle_position;
	};


	class Boundary {
	public:
		Boundary(const double thickness, const bool couples_axis, const bool force_wrap, const bool may_change_particle_pos):
			topology(thickness, couples_axis, force_wrap, may_change_particle_pos)
		{}

		// TODO make this bindable to R-Values
		template<env::FieldMask IncomingMask, env::IsUserData U>
		void invoke_apply(this const auto & self,env::ParticleRef<IncomingMask, U> & particle, const env::Box & domain_box, Face face) noexcept {
			static_assert(
			   requires { { self.apply(particle, domain_box, face) } -> std::same_as<void>; },
			   "BoundaryCondition subclass must implement: void dispatch_apply(particle)"
			);

			using Derived = std::remove_cvref_t<decltype(self)>;

			// check for fields requirements
			static_assert(
				requires { Derived::fields; },
				"Force subclass must define 'static constexpr env::FieldMask fields'"
			);

			constexpr env::FieldMask Required = Derived::fields;

			static_assert(
			  (IncomingMask & Required) == Required,
			  "ParticleView is missing required fields for this Force."
		  );

			self.apply(particle, domain_box, face);
		}

		const Topology topology;
	};



	// define boundary concept
	template <class BC>
	concept IsBoundary = std::derived_from<BC, Boundary>;

	// define boundary pack
	template<IsBoundary... BCs>
	struct BoundaryPack {
	};

	// constrained variable template
	template<class... BCs>
	requires (IsBoundary<BCs> && ...)
	inline constexpr BoundaryPack<BCs...> boundaries {};


	// Concept to check if a type T is a ControllerPack
	template<typename T>
	inline constexpr bool is_boundary_pack_v = false; // Default

	template<IsBoundary... BCs>
	inline constexpr bool is_boundary_pack_v<BoundaryPack<BCs...>> = true; // Specialization

	template<typename T>
	concept IsBoundaryPack = is_boundary_pack_v<std::remove_cvref_t<T>>;


	namespace internal {

		struct BoundarySentinel : Boundary {
			static constexpr env::FieldMask fields = +env::Field::none;

			BoundarySentinel(): Boundary(-1, false, false, false) {}

			template<env::FieldMask IncomingMask, env::IsUserData U>
			void apply(env::ParticleRef<IncomingMask, U> &, const env::Box &, const Face) const noexcept {
				AP_ASSERT(false, "apply called on null boundary! this should never happen");
			}
		};

		// define boundary variant concept
		template<typename T>
		struct is_boundary_variant : std::false_type {};

		template<IsBoundary... BCs>
		struct is_boundary_variant<std::variant<BCs...>> : std::true_type {};

		template<typename T>
		concept IsBoundaryVariant = is_boundary_variant<T>::value;

		template<class... BCs>
		struct VariantType {
			static constexpr bool has_absorb = same_as_any<Open, BCs...>;
			using type = std::conditional_t<
				has_absorb,
				std::variant<BoundarySentinel, BCs...>,
				std::variant<BoundarySentinel, Open, BCs...>
			>;
		};

		template<class... BCs>
		using VariantType_t = typename VariantType<BCs...>::type;
	}
}


