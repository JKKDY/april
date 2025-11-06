#pragma once

#include <variant>

#include "april/boundaries/boundary.h"

namespace april::boundary::internal {

	template<class BVariant>
	const Topology& get_topology(const BVariant& v) {
		return std::visit([]<typename BC>(const BC& bc) -> const Topology& {
			return bc.topology;
		}, v);
	}

	template<class BVariant>
	Topology& get_topology(BVariant& v) {
		return std::visit([]<typename BC>(BC& bc) -> Topology& {
			return bc.topology;
		}, v);
	}



	template<IsBoundaryVariant BVariant> class CompiledBoundary{
	public:
		CompiledBoundary(const BVariant & boundary, const env::Box & boundary_region):
			region(boundary_region), topology(get_topology(boundary)), boundary_v(boundary) {
		}

		template<env::FieldMask M, env::IsUserData UserData>
		void apply(env::ParticleRef<M, UserData> & p, const env::Box & domain_box, const Face face) const noexcept {
			std::visit([&](const auto& bc) {
				bc.dispatch_apply(p, domain_box, face);
			}, boundary_v);
		}

		const env::Box region;
		const Topology topology;
	// private:
		const BVariant boundary_v;
	};


	template<IsBoundaryVariant BVariant>
	CompiledBoundary<BVariant> compile_boundary(const BVariant & boundary, const env::Box & simulation_box, const Face face) {
		constexpr double NEG_INF = std::numeric_limits<double>::lowest() / 4; // divide by 4 to avoid overflow
		constexpr double POS_INF = std::numeric_limits<double>::max() / 4;

		const vec3 POS_INF_VEC = {POS_INF, POS_INF, POS_INF};
		const vec3 NEG_INF_VEC = {NEG_INF, NEG_INF, NEG_INF};

		const int  ax   = axis_of_face(face); // 0:x, 1:y, 2:z (see vec3)
		const bool plus = face_sign_pos(face);

		const Topology& topology = get_topology(boundary);
		const double thickness = topology.boundary_thickness;

		vec3 min, max;

		if (thickness >= 0.0) { // inside the simulation domain
			min = simulation_box.min;
			max = simulation_box.max;

			const double d =  std::clamp(thickness, 0.0, simulation_box.extent[ax]);
			if (plus) {
				min[ax] += max[ax] - d;
			} else {
				max[ax] = min[ax] + d;
			}
		} else { // outside
			if (plus) {
				max = POS_INF_VEC;
				min = NEG_INF_VEC;
				min[ax] = simulation_box.max[ax];
			} else {
				min = NEG_INF_VEC;
				max = POS_INF_VEC / 2;
				max[ax] = simulation_box.min[ax];
			}
		}

		return internal::CompiledBoundary<BVariant>(boundary, env::Box{min, max});
	}



	template<IsBoundaryVariant BVariant>
	struct BoundaryTable {

		BoundaryTable(const std::array<BVariant, 6> & boundaries, const env::Box & simulation_box):
			table({
				compile_boundary<BVariant>(boundaries[face_to_int(Face::XMinus)], simulation_box, Face::XMinus),
				compile_boundary<BVariant>(boundaries[face_to_int(Face::XPlus )], simulation_box, Face::XPlus ),
				compile_boundary<BVariant>(boundaries[face_to_int(Face::YMinus)], simulation_box, Face::YMinus),
				compile_boundary<BVariant>(boundaries[face_to_int(Face::YPlus )], simulation_box, Face::YPlus ),
				compile_boundary<BVariant>(boundaries[face_to_int(Face::ZMinus)], simulation_box, Face::ZMinus),
				compile_boundary<BVariant>(boundaries[face_to_int(Face::ZPlus )], simulation_box, Face::ZPlus ),
			})
		{}

		const CompiledBoundary<BVariant> & get_boundary(const Face face) const {
			return table[face_to_int(face)];
		}

	private:
		std::array<CompiledBoundary<BVariant>, 6> table;
	};
}