#pragma once
#include "core/environment.hpp"

#include "april/boundaries/boundary.hpp"
#include "april/boundaries/absorb.hpp"
#include "april/boundaries/open.hpp"
#include "april/boundaries/periodic.hpp"
#include "april/boundaries/reflective.hpp"
#include "april/boundaries/repulsive.hpp"

#include "april/forces/force.hpp"
#include "april/forces/harmonic.hpp"
#include "april/forces/gravity.hpp"
#include "april/forces/lennard_jones.hpp"
#include "april/forces/no_force.hpp"
#include "april/forces/coulomb.hpp"

#include "april/controllers/controller.hpp"
#include "april/controllers/thermostat.hpp"

#include "april/fields/field.hpp"
#include "april/fields/uniform_field.hpp"
#include "april/fields/local_field.hpp"

#include "april/monitors/monitor.hpp"
#include "april/monitors/terminal_output.hpp"
#include "april/monitors/binary_output.hpp"
#include "april/monitors/progressbar.hpp"
#include "april/monitors/benchmark.hpp"

#include "april/containers/container.hpp"
#include "april/containers/linked_cells/lc_aos.hpp"
#include "april/containers/linked_cells/lc_soa.hpp"
#include "april/containers/linked_cells/lc_aosoa.hpp"

#include "april/containers/linked_cells/lc_config.hpp"
#include "april/containers/cell_orderings.hpp"

#include "april/containers/direct_sum/ds_aos.hpp"
#include "april/containers/direct_sum/ds_soa.hpp"
#include "april/containers/direct_sum/ds_aosoa.hpp"

#include "april/integrators/integrator.hpp"
#include "april/integrators/velocity_verlet.hpp"
#include "april/integrators/yoshida4.hpp"

#include "april/core/build.hpp"
#include "april/core/system.hpp"
#include "april/core/context.hpp"


namespace april {

	// Environment
	// using core::Environment;
	// using ParticleCuboid;
	// using ParticleSphere;
	// using Particle;
	// using core::particle_attributes; // template tag

	// using core::between_types;
	// using core::to_type;
	// using core::between_ids;

	// Boundary
	// using boundary::Face;
	// using boundary::all_faces;

	// using boundary::boundaries; // template pack
	// using boundary::Absorb;
	// using boundary::Open;
	// using boundary::Periodic;
	// using boundary::Reflective;
	// using boundary::Repulsive;

	// Forces
	// using force::forces; // template pack
	// using force::Harmonic;
	// using force::NoForce;
	// using force::Gravity;
	// using force::LennardJones;
	// using force::Coulomb;

	// Controllers
	// using controller::controllers;  // template pack
	// using controller::Controller;
	// using controller::VelocityScalingThermostat;

	// Fields
	// using field::fields; // template pack
	// using field::UniformField;
	// using field::LocalForceField;

	// Containers
	// using container::DirectSumAoS;
	// using container::DirectSumSoA;
	// using container::DirectSumAoSoA;
	// using container::LinkedCellsAoS;
	// using container::LinkedCellsSoA;
	// using container::LinkedCellsAoSoA;

	using container::hilbert_order;
	using container::morton_order;

	// System
	// using core::System;
	// using core::build_system;
	// using core::BuildInfo;

	// Monitors
	// using monitor::monitors;  // template pack
	// using monitor::BinaryOutput;
	// using monitor::TerminalOutput;
	// using monitor::ProgressBar;
	// using monitor::Benchmark;

	// Integrators
	// using integrator::VelocityVerlet;
	// using integrator::Yoshida4;

	// shared
	// using Trigger;
}













