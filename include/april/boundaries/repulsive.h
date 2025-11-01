#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {

	template <typename F, FieldMask M, IsUserData U>
	concept IsBoundaryForce = requires(F f, const ParticleRef<M,U> & p, double dist) {
		{ f.apply(p, dist) } -> std::convertible_to<double>;
		{ f.cutoff() } -> std::convertible_to<double>;
	};



	template <class Force>
	struct Repulsive : Boundary {
		static constexpr FieldMask fields = Field::position | Field::force | Field::old_position;

		explicit Repulsive(Force & force, const bool simulate_halo=false):
		Boundary(force.cutoff(), false, false, false),
		boundary_force(force), simulate_halo(simulate_halo) {}

		template<IsUserData UserData> requires IsBoundaryForce<Force, fields, UserData>
		void apply(ParticleRef<fields, UserData> particle, const Box & domain_box, const Face face) const noexcept{
			const int is_plus = face_sign_pos(face);
			const int ax = axis_of_face(face);

			AP_ASSERT(particle.position[ax] >= domain_box.min[ax] && particle.position[ax] <= domain_box.max[ax],
			"particle should be inside domain on specified axis! \n\t face:"  + std::to_string(face_to_int(face)) +
			"\n\t" + particle.position.to_string() + "  old pos: " + particle.old_position.to_string() );

			const double wall_position  = (is_plus? domain_box.max : domain_box.min)[ax];
			double distance = std::abs(wall_position  - particle.position[ax]);

			// simulate as if the particle where interacting with itself on the other side of the boundary
			if (simulate_halo) distance  *= 2;

			// for repulsive we always assume boundary_force returns a positive scalar
			double magnitude = boundary_force.apply(particle, distance );
			int direction = is_plus ? -1 : 1;  // push inward
			particle.force[ax] += direction * magnitude;
		}
	private:
		Force boundary_force;
		bool simulate_halo;
	};



	struct ExponentialForce {
		double A;       // amplitude
		double lambda;  // decay length
		double rc;      // cutoff distance

		ExponentialForce(const double A, const double lambda, const double rc)
			: A(A), lambda(lambda), rc(rc) {}

		[[nodiscard]] double cutoff() const noexcept { return rc; }

		template <FieldMask M, IsUserData U>
		[[nodiscard]] double apply(const ParticleRef<M, U>, const double distance) const noexcept {
			if (distance > rc) return 0.0;
			return A * std::exp(-distance / lambda);
		}
	};

	struct PowerLawForce {
		double A;
		double n;
		double rc;

		PowerLawForce(const double A, const double n, const double rc)
			: A(A), n(n), rc(rc) {}

		[[nodiscard]] double cutoff() const noexcept { return rc; }

		template <FieldMask M, IsUserData U>
		[[nodiscard]] double apply(const ParticleRef<M, U>, const double distance) const noexcept {
			if (distance > rc) return 0.0;
			return A / std::pow(distance, n);
		}
	};

	struct LennardJones93Force {
		double epsilon;
		double sigma;
		double rc;

		LennardJones93Force(const double epsilon, const double sigma, const double rc)
			: epsilon(epsilon), sigma(sigma), rc(rc) {}

		[[nodiscard]] double cutoff() const noexcept { return rc; }

		template <FieldMask M, IsUserData U>
		[[nodiscard]] double apply(const ParticleRef<M, U>, const double distance) const noexcept {
			if (distance > rc) return 0.0;
			const double sr = sigma / distance;
			const double sr3 = sr * sr * sr;
			const double sr9 = sr3 * sr3 * sr3;
			return 4.0 * epsilon * (3.0 * sr3 - 9.0 * sr9);
		}
	};


	struct AdhesiveLJForce {
		double epsilon;
		double sigma;
		double rc;

		AdhesiveLJForce(const double epsilon, const double sigma, const double rc)
			: epsilon(epsilon), sigma(sigma), rc(rc) {}

		[[nodiscard]] double cutoff() const noexcept { return rc; }

		template <FieldMask M, IsUserData U>
		[[nodiscard]] double apply(const ParticleRef<M, U>, const double distance) const noexcept {
			if (distance > rc) return 0.0;
			const double sr = sigma / distance;
			const double sr6 = std::pow(sr, 6);
			const double sr12 = sr6 * sr6;
			const double f = 24.0 * epsilon * (2.0 * sr12 - sr6) / distance;
			return std::abs(f);  // ensure positive scalar magnitude
		}
	};



}