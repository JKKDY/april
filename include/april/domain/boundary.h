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

	constexpr int to_int(Face f) noexcept {
		return static_cast<int>(f);
	}

	enum class Side: uint8_t {
		Inside=1, Outside=2, Both=3
	};

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

		void dispatch_apply(this auto&& self, env::impl::Particle & particle) noexcept {
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

	// Accepts: std::variant<std::monostate, BCs...>
	template<typename T>
	struct is_boundary_variant : std::false_type {};

	template<typename... BCs>
	struct is_boundary_variant<std::variant<std::monostate, BCs...>>
	: std::bool_constant<(IsBoundary<BCs> && ...)> {};

	// define boundary variant concept
	template<typename T>
	concept IsBoundaryVariant = is_boundary_variant<T>::value;

	// define boundary pack
	template<IsBoundary... BCs>
	struct BoundaryPack {};

	template<class... BCs>
	requires (IsBoundary<BCs> && ...)
	inline constexpr BoundaryPack<BCs...> boundaries {};





	namespace impl {

		template<class BVariant>
		const Topology& get_topology(const BVariant& v) {
			return std::visit([](const auto& bc) -> const Topology& {
				using T = std::decay_t<decltype(bc)>;
				if constexpr (requires(const T& x) { x.topology; }) {
					return bc.topology;
				} else {
					throw std::logic_error("Method called on Boundary Variant in std::monostate");
				}
			}, v);
		}

		template<class BVariant>
		Topology& get_topology(BVariant& v) {
			return std::visit([](auto& bc) -> Topology& {
				using T = std::decay_t<decltype(bc)>;
				if constexpr (requires(T& x) { x.topology; }) {
					return bc.topology;
				} else {
					throw std::logic_error("Method called on Boundary Variant in std::monostate");
				}
			}, v);
		}


		template<IsBoundaryVariant BVariant> class CompiledBoundary{
		public:
			CompiledBoundary(const BVariant & boundary, const Domain & boundary_region):
				region(boundary_region), boundary_v(boundary) {

				std::visit([this](const auto & bc) {
					using T = std::decay_t<decltype(bc)>; // get the type of the alternative
					if constexpr (requires(T& x) { x.dispatch_apply; }) {
						apply_fn = &thunk<T>; // store the thunk in apply_fn
					}
				}, boundary_v);
			}

			void apply(Particle & p) const noexcept {
				apply_fn(this, p);
			}

			const Domain region;
		private:
			template<typename T>
			static void thunk(const CompiledBoundary * self, Particle & p ) noexcept {
				// thunk for uniform function pointer type regardless of underlying variant type
				// get the alternative from the variant and call apply on the particle
				std::get<T>(self->boundary_v).dispatch_apply(p);
			}

			using ApplyFn = void (*)(const CompiledBoundary*, Particle&);
			ApplyFn apply_fn = nullptr;

			const BVariant boundary_v;
		};

		template<IsBoundaryVariant BVariant>
		CompiledBoundary<BVariant> compile_boundary(const BVariant & boundary, const Domain & env_domain, Face face) {

			constexpr double NEG_INF = std::numeric_limits<double>::lowest();
			constexpr double POS_INF = std::numeric_limits<double>::max();

			const auto axis_of  = [](Face f) noexcept -> int  { return static_cast<int>(f) / 2; };
			const auto is_plus  = [](Face f) noexcept -> bool { return (static_cast<int>(f) & 1) != 0; };
			const auto in_width = [](const double t, const double L) noexcept -> double {
				// clamp to [0, L]
				const double d = std::max(0.0, t);
				return std::min(d, L);
			};

			Domain region = env_domain;  // start with full domain, override later
			const int  ax  = axis_of(face); // 0:x, 1:y, 2:z (see vec3)
			const bool plus = is_plus(face);

			const Topology& topo = get_topology(boundary);
			const double t = topo.boundary_thickness;

			if (t >= 0.0) {
				const double d = in_width(t, env_domain.extent[ax]);
				region.extent[ax] = d;
				region.origin[ax] = plus
					? (env_domain.origin[ax] + (env_domain.extent[ax] - d))   // [max-d, max]
					:  env_domain.origin[ax];                                 // [min, min+d]
			} else {
				if (plus) {
					const double edge = env_domain.origin[ax] + env_domain.extent[ax];
					region.origin[ax] = edge;                 // [edge, +MAX]
					region.extent[ax] = POS_INF - edge;       // so origin+extent == +MAX
				} else {
					region.origin[ax] = NEG_INF;              // [LOWEST, edge]
					region.extent[ax] = env_domain.origin[ax] - NEG_INF; // so sum == edge
				}
			}

			return impl::CompiledBoundary<BVariant>(boundary, region);
		}

		template<IsBoundaryVariant BVariant>
		struct BoundaryTable {

			BoundaryTable(const std::array<BVariant, 6> & boundaries_, const Domain & env_domain):
				boundaries({
					compile_boundary<BVariant>(boundaries_[to_int(Face::XMinus)], env_domain, Face::XMinus),
					compile_boundary<BVariant>(boundaries_[to_int(Face::XPlus)], env_domain, Face::XPlus),
					compile_boundary<BVariant>(boundaries_[to_int(Face::YMinus)], env_domain, Face::YMinus),
					compile_boundary<BVariant>(boundaries_[to_int(Face::YPlus)], env_domain, Face::YPlus),
					compile_boundary<BVariant>(boundaries_[to_int(Face::ZMinus)], env_domain, Face::ZMinus),
					compile_boundary<BVariant>(boundaries_[to_int(Face::ZPlus)], env_domain, Face::ZPlus),
				})
			{}

			CompiledBoundary<BVariant> & boundary(const Face face) {
				return boundaries[to_int(face)];
			}

			std::array<CompiledBoundary<BVariant>, 6> boundaries;
		};

	}



	struct Absorb : Boundary {
		Absorb(): Boundary(-1, false, false) {}

		void apply(impl::Particle & particle) const noexcept{
			particle.state = ParticleState::DEAD;
		}
	};

	struct Outflow : Boundary {
		Outflow(): Boundary(-1, false, false) {}

	};

	struct Periodic : Boundary {
		Periodic(): Boundary(-1, true, true) {}

	};

	struct Reflective : Boundary {
		Reflective(): Boundary(-1, false, false) {}
	};

	template <class Force>
	struct Repulsive : Boundary {
		explicit Repulsive(Force & force): Boundary(force.cutoff(), false, false) {}

	};
}


