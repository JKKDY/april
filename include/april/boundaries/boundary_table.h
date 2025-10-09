#pragma once

#include <variant>

#include "april/boundaries/boundary.h"

namespace april::boundary::internal {

		template<class BVariant>
		const Topology& get_topology(const BVariant& v) {
			return std::visit([]<typename BC>(const BC& bc) -> const Topology& {
				using T = std::decay_t<BC>;
				if constexpr (requires(const T& x) { x.topology; }) {
					return bc.topology;
				} else {
					throw std::logic_error("Method called on Boundary Variant in std::monostate");
				}
			}, v);
		}

		template<class BVariant>
		Topology& get_topology(BVariant& v) {
			return std::visit([]<typename BC>(BC& bc) -> Topology& {
				using T = std::decay_t<BC>;
				if constexpr (requires(T& x) { x.topology; }) {
					return bc.topology;
				} else {
					throw std::logic_error("Method called on Boundary Variant in std::monostate");
				}
			}, v);
		}


		template<BoundaryVariant BVariant> class CompiledBoundary{
		public:
			CompiledBoundary(const BVariant & boundary, const env::Domain & boundary_region):
				region(boundary_region), boundary_v(boundary) {

				std::visit([this]<typename BC>(const BC &) {
					using T = std::decay_t<BC>; // get the type of the alternative
					if constexpr (requires(T& x) { x.dispatch_apply; }) {
						apply_fn = &thunk<T>; // store the thunk in apply_fn
					}
				}, boundary_v);
			}

			void apply(env::internal::Particle & p) const noexcept {
				apply_fn(this, p);
			}

			const env::Domain region;
		private:
			template<typename T>
			static void thunk(const CompiledBoundary * self, env::internal::Particle & p ) noexcept {
				// thunk for uniform function pointer type regardless of underlying variant type
				// get the alternative from the variant and call apply on the particle
				std::get<T>(self->boundary_v).dispatch_apply(p);
			}

			using ApplyFn = void (*)(const CompiledBoundary*, env::internal::Particle&);
			ApplyFn apply_fn = nullptr;

			const BVariant boundary_v;
		};



		template<BoundaryVariant BVariant>
		CompiledBoundary<BVariant> compile_boundary(const BVariant & boundary, const env::Domain & env_domain, const Face face) {

			constexpr double NEG_INF = std::numeric_limits<double>::lowest();
			constexpr double POS_INF = std::numeric_limits<double>::max();

			const auto axis_of  = [](Face f) noexcept -> int  { return static_cast<int>(f) / 2; };
			const auto is_plus  = [](Face f) noexcept -> bool { return (static_cast<int>(f) & 1) != 0; };
			const auto in_width = [](const double t, const double L) noexcept -> double {
				// clamp to [0, L]
				const double d = std::max(0.0, t);
				return std::min(d, L);
			};

			env::Domain region = env_domain;  // start with full domain, override later
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

			return internal::CompiledBoundary<BVariant>(boundary, region);
		}



		template<BoundaryVariant BVariant>
		struct BoundaryTable {

			BoundaryTable(const std::array<BVariant, 6> & boundaries, const env::Domain & env_domain):
				table({
					compile_boundary<BVariant>(boundaries[to_int(Face::XMinus)], env_domain, Face::XMinus),
					compile_boundary<BVariant>(boundaries[to_int(Face::XPlus )], env_domain, Face::XPlus ),
					compile_boundary<BVariant>(boundaries[to_int(Face::YMinus)], env_domain, Face::YMinus),
					compile_boundary<BVariant>(boundaries[to_int(Face::YPlus )], env_domain, Face::YPlus ),
					compile_boundary<BVariant>(boundaries[to_int(Face::ZMinus)], env_domain, Face::ZMinus),
					compile_boundary<BVariant>(boundaries[to_int(Face::ZPlus )], env_domain, Face::ZPlus ),
				})
			{}

			CompiledBoundary<BVariant> & get_boundary(const Face face) {
				return table[to_int(face)];
			}

		private:
			std::array<CompiledBoundary<BVariant>, 6> table;
		};
}