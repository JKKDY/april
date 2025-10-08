#pragma once
#include "april/env/environment.h"
#include "april/env/particle.h"

#include "april/boundaries/boundary.h"
#include "april/boundaries/absorb.h"
#include "april/boundaries/outflow.h"
#include "april/boundaries/periodic.h"
#include "april/boundaries/reflective.h"
#include "april/boundaries/repulsive.h"

#include "april/forces/force.h"
#include "april/forces/harmonic.h"
#include "april/forces/inverse_square.h"
#include "april/forces/lennard_jones.h"
#include "april/forces/no_force.h"

#include "april/monitors/monitor.h"
#include "april/monitors/binary_output.h"
#include "april/monitors/progressbar.h"
#include "april/monitors/benchmark.h"

#include "april/containers/container.h"
#include "april/containers/direct_sum.h"
#include "april/containers/linked_cells.h"

#include "april/integrators/integrator.h"
#include "april/integrators/stoermer_verlet.h"
#include "april/integrators/yoshida4.h"

#include "april/boundaries/boundary.h"



#define APRIL_API

// TODO rework namespace
// TODO use modules
namespace april {

	// Environment
	using env::Environment;

	using env::ParticleCuboid;
	using env::ParticleSphere;
	using env::ParticleState;
	using env::Particle;
	using env::ParticleID;
	using env::ParticleType;
	using env::PARTICLE_ID_DONT_CARE;

	using env::between_types;
	using env::to_type;
	using env::between_ids;

	// Boundary
	using boundary::boundaries;
	using boundary::Absorb;
	using boundary::Outflow;
	using boundary::Periodic;
	using boundary::Reflective;
	using boundary::Repulsive;

	// Forces
	using force::forces;
	using force::Harmonic;
	using force::NoForce;
	using force::InverseSquare;
	using force::LennardJones;

	// Containers
	using container::DirectSum;
	using container::LinkedCells;

	
	// System
	using core::System;
	using core::build_system;
	using core::UserToInternalMappings;

	// Monitors
	using monitor::monitors;
	using monitor::BinaryOutput;
	using monitor::TerminalOutput;
	using monitor::ProgressBar;
	using monitor::Benchmark;

	// Integrators
	using integrator::StoermerVerlet;
	using integrator::Yoshida4;

	using env::impl::ParticleView; // move to non impl


	namespace ext {
		using monitor::Monitor;
		using boundary::Boundary;
		using container::impl::Container;
		using integrator::impl::Integrator;
		using env::impl::Particle;
	}
}