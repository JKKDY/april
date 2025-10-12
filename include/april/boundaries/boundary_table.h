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
				region(boundary_region), topology(get_topology(boundary)), boundary_v(boundary) {

				std::visit([&](auto const& bc) {
					using T = std::decay_t<decltype(bc)>; // get the type of the alternative

					if constexpr (std::same_as<T, std::monostate>) {
						throw std::logic_error("Trying to set the boundary thunk on a monostate should never happen!");
					} else {
						static_assert(IsBoundary<T>, "Variant alternative must satisfy IsBoundary");
						apply_fn =&thunk<T>; // store the thunk in apply_fn
					}
				}, boundary_v);
			}

			void apply(env::internal::Particle & p, const env::Box & domain_box, Face face) const noexcept {
				apply_fn(this, p, domain_box, face);
			}

			const env::Domain region;
			const Topology topology;
		private:
			template<typename T>
			static void thunk(const CompiledBoundary * self, env::internal::Particle & p, const env::Box & domain_box, Face face) noexcept {
				// thunk for uniform function pointer type regardless of underlying variant type
				// get the alternative from the variant and call apply on the particle
				std::get<T>(self->boundary_v).dispatch_apply(p, domain_box, face);
			}

			using ApplyFn = void (*)(const CompiledBoundary*, env::internal::Particle&, const env::Box&, Face);
			ApplyFn apply_fn = nullptr;

			const BVariant boundary_v;
		};


		template<BoundaryVariant BVariant>
		CompiledBoundary<BVariant> compile_boundary(const BVariant & boundary, const env::Domain & env_domain, const Face face) {
			// TODO we expect env_domain.extent to be positive on all axis -> assert this

			constexpr double NEG_INF = std::numeric_limits<double>::lowest() / 2; // divide by half to avoid overflow
			constexpr double POS_INF = std::numeric_limits<double>::max() / 2;

			const vec3 POS_INF_VEC = {POS_INF, POS_INF, POS_INF};
			const vec3 NEG_INF_VEC = {NEG_INF, NEG_INF, NEG_INF};

			const int  ax   = axis_of_face(face); // 0:x, 1:y, 2:z (see vec3)
			const bool plus = face_sign_pos(face);

			const Topology& topo = get_topology(boundary);
			const double thickness = topo.boundary_thickness;

			env::Domain region = {};

			if (thickness >= 0.0) { // inside the simulation domain
				region.extent = env_domain.extent;
				region.origin = env_domain.origin;

				const double d =  std::clamp(thickness, 0.0, env_domain.extent[ax]);
				if (plus) region.origin[ax] += env_domain.extent[ax] - d;
				region.extent[ax] = d;
			} else { // outside
				if (plus) {
					region.extent = POS_INF_VEC;
					region.origin = NEG_INF_VEC / 2;
					region.origin[ax] = env_domain.origin[ax] + env_domain.extent[ax];
				} else {
					region.extent = NEG_INF_VEC;
					region.origin = POS_INF_VEC / 2;
					region.origin[ax] = env_domain.origin[ax];
				}
			}

			return internal::CompiledBoundary<BVariant>(boundary, region);
		}



		template<BoundaryVariant BVariant>
		struct BoundaryTable {

			BoundaryTable(const std::array<BVariant, 6> & boundaries, const env::Domain & env_domain):
				table({
					compile_boundary<BVariant>(boundaries[face_to_int(Face::XMinus)], env_domain, Face::XMinus),
					compile_boundary<BVariant>(boundaries[face_to_int(Face::XPlus )], env_domain, Face::XPlus ),
					compile_boundary<BVariant>(boundaries[face_to_int(Face::YMinus)], env_domain, Face::YMinus),
					compile_boundary<BVariant>(boundaries[face_to_int(Face::YPlus )], env_domain, Face::YPlus ),
					compile_boundary<BVariant>(boundaries[face_to_int(Face::ZMinus)], env_domain, Face::ZMinus),
					compile_boundary<BVariant>(boundaries[face_to_int(Face::ZPlus )], env_domain, Face::ZPlus ),
				})
			{}

			CompiledBoundary<BVariant> & get_boundary(const Face face) {
				return table[face_to_int(face)];
			}

		private:
			std::array<CompiledBoundary<BVariant>, 6> table;
		};
}