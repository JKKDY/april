#pragma once

#include <cstdint>
#include "april/common.h"
#include "april/domain/domain.h"

namespace april::env {

	enum class Face : uint8_t {
		XMinus = 0, XPlus = 1,
		YMinus = 2, YPlus = 3,
		ZMinus = 4, ZPlus = 5,
	};

	enum class Side: uint8_t {
		Inside=1, Outside=2, Both=3
	};


	struct Topology {
		Side side;
		double thickness; // negative means the implementation decide
		bool couples_axis;
		bool force_wrap = false;
	};

	class BoundaryCondition {
	public:
		BoundaryCondition(const Side side, const double thickness, const bool couples_axis, const bool force_wrap):
			topology_(side, thickness, couples_axis, force_wrap)
		{}

		void dispatch_apply(this auto&& self, env::impl::Particle & particle) noexcept {
			static_assert(
			   requires { { self.apply(particle) } -> std::same_as<void>; },
			   "BoundaryCondition subclass must implement: void dispatch_apply(particle)"
			);

			self.apply(particle);
		}

		[[nodiscard]] Topology topology() const {
			return topology_;
		}

	private:
		Topology topology_;
	};


	// boundary concept
	template <class BC>
	concept IsBoundary = std::derived_from<BC, BoundaryCondition>;


	// boundary variant
	template<typename T>
	struct is_boundary_variant : std::false_type {};

	template<IsBoundary... BCs>
    struct is_boundary_variant<std::variant<BCs...>> : std::true_type {};

	template<typename T>
	concept BoundaryVariant = is_boundary_variant<T>::value;


	// boundary pack
	template<IsBoundary... BCs>
	struct BoundaryPack { /* empty */ };

	template<class... BCs>
	inline constexpr BoundaryPack<BCs...> boundaries {};





	template <IsBoundary... BCs>
	struct BoundaryImpl {

		void apply(env::impl::Particle & particle) const noexcept {
			boundary_condition.dispatch_apply(particle);
		}


	private:
		std::variant<BCs...> boundary_condition;
		Domain region_;
	};

	template<IsBoundary... BCs>
	struct BoundaryTable {
		using boundary_variant_t = std::variant<BoundaryImpl<BCs...>>;
		// std::array<boundary_variant_t, 6>;
	};


	class Absorb : public BoundaryCondition {

	};

	class Outflow : public BoundaryCondition {

	};

	class Periodic : public BoundaryCondition {

	};

	class Reflective : public BoundaryCondition {

	};

	class Repulsive : public BoundaryCondition {

	};

	// class Absorb : public BoundaryCondition{
	// public:
	// 	Absorb() {
	//
	// 	}
	//
	// 	void apply(env::impl::Particle & particle) const noexcept{
	// 		particle.state = env::ParticleState::DEAD;
	// 	}
	// };
	//
	// class Periodic : public BoundaryCondition{
	// public:
	// 	void apply(env::impl::Particle & particle) const noexcept{
	// 		// particle.position =
	// 	}
	// private:
	// 	env::Domain sim_domain;
	// };
}


