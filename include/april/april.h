#pragma once
#include "april/env/environment.h"
#include "april/env/particle.h"
#include "april/env/interaction.h"

#include "april/io/monitor.h"
#include "april/io/output.h"

#include "april/containers/container.h"
#include "april/containers/direct_sum.h"
#include "april/containers/linked_cells.h"


#include "april/core/integrator.h"
#include "april/core/stoermer_verlet.h"


// Define APRIL_API for DLL import/export
// #ifdef _WIN32
//   #ifdef APRIL_BUILD_DLL
//     #define APRIL_API __declspec(dllexport)
//   #else
//     #define APRIL_API __declspec(dllimport)
//   #endif
// #else
//   #define APRIL_API
// #endif

#define APRIL_API


namespace april {

	using core::DirectSum;
	using core::LinkedCells;

	using io::BinaryOutput;
	using io::TerminalOutput;

	using env::Environment;

	using env::ParticleCuboid;
	using env::ParticleSphere;
	using env::ParticleState;
	using env::Particle;
	using env::ParticleID;
	using env::ParticleType;
	using env::PARTICLE_ID_UNDEFINED;

	using env::Force;
	using env::Harmonic;
	using env::NoForce;
	using env::InverseSquare;
	using env::LennardJones;

	using utils::Vec3;

	template <io::IsMonitor... TMonitors> using StoermerVerlet = core::StoermerVerlet<TMonitors...>;

	namespace ext {
		using io::Monitor;
		using core::Container;
		template <io::IsMonitor... TMonitors> using Integrator = core::impl::Integrator<TMonitors...>;
	}
}