#pragma once
#include "april/env/environment.h"
#include "april/env/particle.h"
#include "april/env/interaction.h"

#include "april/io/monitor.h"
#include "april/io/output.h"
#include "april/io/status.h"
#include "april/io/performance.h"

#include "april/containers/container.h"
#include "april/containers/direct_sum.h"
#include "april/containers/linked_cells.h"

#include "april/core/integrator.h"
#include "april/core/stoermer_verlet.h"
#include "april/core/yoshida4.h"



#define APRIL_API

// TODO rework namespace
// TODO use modules
namespace april {
	using env::Environment;

	using env::ParticleCuboid;
	using env::ParticleSphere;
	using env::ParticleState;
	using env::Particle;
	using env::ParticleID;
	using env::ParticleType;
	using env::PARTICLE_ID_DONT_CARE;

	using env::Harmonic;
	using env::NoForce;
	using env::InverseSquare;
	using env::LennardJones;
	using env::forces;

	using env::between_types;
	using env::to_type;
	using env::between_ids;

	using cont::DirectSum;
	using cont::LinkedCells;

	using core::System;
	using core::build_system;
	using core::UserToInternalMappings;

	using io::BinaryOutput;
	using io::TerminalOutput;
	using io::ProgressBar;
	using io::Benchmark;
	using io::monitors;

	using core::StoermerVerlet;
	using core::Yoshida4;

	namespace ext {
		using io::Monitor;
		using cont::impl::Container;
		using core::impl::Integrator;
	}
}