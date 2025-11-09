#pragma once
#include "april/env/environment.h"
#include "april/particle/particle_fields.h"

#include "april/boundaries/boundary.h"
#include "april/boundaries/absorb.h"
#include "april/boundaries/open.h"
#include "april/boundaries/periodic.h"
#include "april/boundaries/reflective.h"
#include "april/boundaries/repulsive.h"

#include "april/forces/force.h"
#include "april/forces/harmonic.h"
#include "april/forces/power_law.h"
#include "april/forces/lennard_jones.h"
#include "april/forces/no_force.h"

#include "april/controllers/controller.h"
#include "april/controllers/thermostat.h"

#include "april/fields/field.h"
#include "april/fields/uniform_field.h"
#include "april/fields/local_field.h"

#include "april/monitors/monitor.h"
#include "april/monitors/terminal_output.h"
#include "april/monitors/binary_output.h"
#include "april/monitors/progressbar.h"
#include "april/monitors/benchmark.h"

#include "april/containers/container.h"
#include "april/containers/direct_sum.h"
#include "april/containers/linked_cells.h"

#include "april/integrators/integrator.h"
#include "april/integrators/stoermer_verlet.h"
#include "april/integrators/yoshida4.h"

#include "april/system/build.h"
#include "april/system/system.h"
#include "april/system/context.h"


#define APRIL_API

namespace april {

	// Environment
	using env::Environment;

	using env::ParticleCuboid;
	using env::ParticleSphere;
	using env::ParticleState;
	using env::ParticleView;
	using env::Particle;
	using env::ParticleID;
	using env::ParticleType;
	using env::particle_data;

	using env::between_types;
	using env::to_type;
	using env::between_ids;

	// Boundary
	using boundary::Face;
	using boundary::all_faces;

	using boundary::boundaries;
	using boundary::Boundary;
	using boundary::Absorb;
	using boundary::Open;
	using boundary::Periodic;
	using boundary::Reflective;
	using boundary::Repulsive;

	// Forces
	using force::forces;
	using force::Harmonic;
	using force::NoForce;
	using force::PowerLaw;
	using force::LennardJones;

	// Controllers
	using controller::controllers;
	using controller::Controller;
	using controller::VelocityScalingThermostat;

	// Fields
	using field::fields;
	using field::UniformField;
	using field::LocalForceField;

	// Containers
	using container::Container;
	using container::DirectSum;
	using container::LinkedCells;

	// System
	using core::System;
	using core::build_system;
	using core::BuildInfo;
	using core::SystemContext;

	// Monitors
	using monitor::monitors;
	using monitor::Monitor;
	using monitor::BinaryOutput;
	using monitor::TerminalOutput;
	using monitor::ProgressBar;
	using monitor::Benchmark;

	// Integrators
	using integrator::Integrator;
	using integrator::StoermerVerlet;
	using integrator::Yoshida4;

	// shared
	using shared::Trigger;
}