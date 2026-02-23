#pragma once

#include "april/boundaries/boundary.hpp"
#include "april/forces/force.hpp"


namespace april::core {

	//--------------
	// EXECUTE BATCH
	//--------------
	template <class ContainerDecl, core::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<
		ContainerDecl, Traits>
	template <ParticleField M, ParallelPolicy P, VectorPolicy V, container::IsBatch Batch, exec::IsKernel Kernel>
	void System<ContainerDecl, Traits>::execute_batch_kernel(const Batch& batch, Kernel&& kernel)  {
			using namespace april::exec::internal;
			constexpr VectorTrait batch_capabilities  = Batch::vector_trait;
			constexpr ExecutionMode kernel_capabilities = Kernel::Mode;

			// map VectorPolicy -> VectorTrait (execution mode)
			constexpr ExecutionMode exec_mode = [] {
				if constexpr (V == VectorPolicy::Vector) {
					return ExecutionMode::Vector; // Force Vector
				}
				else if constexpr (V == VectorPolicy::Scalar) {
					return ExecutionMode::Scalar; // Force Scalar
				}
				else { // Auto
					return kernel_capabilities; // return kernel capabilities. Batch will figure out what to do
				}
			}();

			// Enforce Contracts
			if constexpr (V == VectorPolicy::Vector) {  // Case user demands strict vectorization
				// Batch must have a vector only path
				static_assert(has_flag(batch_capabilities, VectorTrait::VectorOnly),
					"Policy Violation: VectorPolicy::Vector requires a fully vectorized Batch. ");

				// Kernel must be Vector-Capable
				static_assert(has_flag(kernel_capabilities, ExecutionMode::Vector),
					"Policy Violation: VectorPolicy::Vector requires a SIMD-capable Kernel.");
			}
			else if constexpr (V == VectorPolicy::Scalar) {  // Case user demands strict scalar computation
				// Batch must have a scalar only path
				static_assert(has_flag(batch_capabilities, VectorTrait::ScalarOnly),
					"Policy Violation: VectorPolicy::Scalar requires a fully scalarized Batch. ");

				//  Kernel must be Scalar-Capable
				static_assert(has_flag(kernel_capabilities, ExecutionMode::Scalar),
					"Policy Violation: VectorPolicy::Scalar requires a Scalar-capable Kernel.");
			}
			else if constexpr (V == VectorPolicy::Auto) {  // Case auto (best effort vectorization)
				constexpr bool CanRunVector = has_flag(batch_capabilities, VectorTrait::VectorOnly) &&
					has_flag(kernel_capabilities, ExecutionMode::Vector);

				constexpr bool CanRunScalar = has_flag(batch_capabilities, VectorTrait::ScalarOnly) &&
					has_flag(kernel_capabilities, ExecutionMode::Scalar);

				constexpr bool CanRunHybrid = has_flag(batch_capabilities, VectorTrait::Mixed) &&
					kernel_capabilities == (ExecutionMode::Vector | ExecutionMode::Scalar);

				static_assert(CanRunVector || CanRunScalar || CanRunHybrid,
					"Compatibility Failure: No valid execution path found between Batch and Kernel capability sets.");
			}

			// Execute
			if constexpr (container::IsBatchAtom<Batch>) {
				batch.template for_each_pair<M, P, exec_mode>(kernel);
			}
			else if constexpr (container::IsBatchAtomRange<Batch>) {
				for (const auto& atom : batch) {
					atom.template for_each_pair<M, P, exec_mode>(kernel);
				}
			}
		}



	//--------------
	// UPDATE FORCES
	//--------------
	template <class C, core::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::update_forces() {

		// batch update lambda. passed into container::for_each_interaction_batch
		auto update_batch = [&]<container::IsBatch Batch, container::IsBCP BCP>(const Batch& batch, BCP && apply_bcp) {

			auto apply_batch_update =  [&] <force::IsForce ForceT> (const ForceT & force) {
				constexpr ParticleField M = ForceT::fields | ParticleField::force | ParticleField::position;

				auto kernel = [&]<bool is_packed>(auto & p1, auto & p2) {
					auto diff = p2.position - p1.position;

					const auto r = [&] {
						if constexpr (std::is_same_v<std::decay_t<BCP>, container::NoBatchBCP>) {
							return diff;
						} else {
							return apply_bcp(diff);
						}
					}();

					if constexpr (is_packed) {
						auto mask = r.norm_squared() > force.cutoff2();
						if (all(mask)) return;
						auto f = force(p1, p2, r);
						auto f_masked = pvec3 {
							select(mask, packed(0), f.x),
							select(mask, packed(0), f.y),
							select(mask, packed(0), f.z)
						};
						p1.force += f_masked;
						p2.force -= f_masked;
					} else {
						if (r.norm_squared() > force.cutoff2()) {
							return;
						}
						vec3 f = force(p1.to_view(), p2.to_view(), r);
						p1.force += f;
						p2.force -= f;
					}
				};

				execute_batch_kernel<M, ParallelPolicy::Serial, VectorPolicy::Auto>(batch, april::universal_kernel(kernel));
			};

			auto [t1, t2] = batch.types;
			force_table.dispatch(t1, t2, apply_batch_update);
		};

		particle_container.template for_each_particle<ParticleField::force>(april::universal_kernel(
			[](auto && p) { p.force = {}; } // reset forces
		));
		particle_container.invoke_for_each_interaction_batch(update_batch);




		// handle id interactions
		auto update_global_batch = [&](const auto & batch) {

			auto apply_batch_update = [&] <force::IsForce ForceT> (const ForceT & force) {
				constexpr ParticleField M = ForceT::fields | ParticleField::force | ParticleField::position;

				for (const auto & [id1, id2] : batch.pairs) {
					auto && p1 = particle_container.template restricted_at_id<M>(id1);
					auto && p2 = particle_container.template restricted_at_id<M>(id2);

					vec3 r = p2.position - p1.position;
					const vec3 f = force(p1.to_view(), p2.to_view(), r);

					p1.force += f;
					p2.force -= f;
				}
			};

			force_table.dispatch_id(batch.id1, batch.id2, apply_batch_update);
		};

		particle_container.invoke_for_each_topology_batch(update_global_batch);
	}



	//-----------------
	// APPLY BOUNDARIES
	//-----------------
	template <class C, core::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::apply_boundary_conditions() {
		particles_to_update_buffer.clear();

		const core::Box domain_box = this->box();

		for (boundary::Face face : boundary::all_faces) {

			const auto& compiled_boundary = boundary_table[face];

			std::vector<size_t> particle_ids = particle_container.invoke_collect_indices_in_region(compiled_boundary.boundary_region);

			auto boundary_condition_inside = [&]<typename B>(const B & bc) {
				constexpr ParticleField M = std::decay_t<B>::fields;

				for (auto p_idx : particle_ids) {
					auto p = particle_container.template at<M>(p_idx);
					bc.apply(p, domain_box, face);

					if (compiled_boundary.topology.may_change_particle_position) {
						particles_to_update_buffer.push_back(p_idx);
					}
				}
			};

			auto boundary_condition_outside = [&]<typename B>(const B & bc) {
				static constexpr ParticleField detect_mask = ParticleField::position | ParticleField::old_position;
				constexpr ParticleField M = std::decay_t<B>::fields | detect_mask;

				for (auto p_idx : particle_ids) {
					auto particle = particle_container.template at<M>(p_idx);

					// make sure the particle exited through the current boundary face
					// solve for intersection of the particles path with the boundary face
					// with the equation y = t * diff + p where:
					// diff is the path traveled, p is the particles starting position and y is the face

					const int ax = axis_of_face(face);
					const vec3 diff = particle.position - particle.old_position;
					const double y = diff[ax] < 0 ? domain_box.min[ax] : domain_box.max[ax];
					const double t = (y - particle.old_position[ax]) / diff[ax];

					const vec3 intersection = t * diff + particle.old_position;

					// and check if that point is on the domains surface
					auto [ax1, ax2] = non_face_axis(face);
					if (domain_box.max[ax1] >= intersection[ax1] && domain_box.min[ax1] <= intersection[ax1] &&
						domain_box.max[ax2] >= intersection[ax2] && domain_box.min[ax2] <= intersection[ax2]) {

						bc.apply(particle, domain_box, face);

						if (compiled_boundary.topology.may_change_particle_position) {
							particles_to_update_buffer.push_back(p_idx);
						}
						}
				}
			};

			if (compiled_boundary.topology.boundary_thickness >= 0) {
				compiled_boundary.dispatch(boundary_condition_inside);
			} else {
				compiled_boundary.dispatch(boundary_condition_outside);
			}
		}

		if (!particles_to_update_buffer.empty()) {
			particle_container.invoke_notify_moved(particles_to_update_buffer);
		}
	}


	//------------------
	// APPLY CONTROLLERS
	//------------------
	template <class C, core::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::apply_controllers() {
		controllers.for_each_item([this](auto & controller) {
			if (controller.should_trigger(trig_context)) {
				controller.dispatch_apply(system_context);
			}
		});
	}


	//-------------
	// APPLY FIELDS
	//-------------
	template <class C, core::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::apply_force_fields() {
		fields.for_each_item([this]<typename F>(F & field) {
			for (size_t i = 0; i < size(); ++i) {
				constexpr ParticleField M = F::fields;
				auto restricted = particle_container.template restricted_at<M>(i);
				field.dispatch_apply(restricted);
			}
		});
	}


	//-------
	// UPDATE
	//-------
	template <class C, core::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::update_all_components() {
		fields.for_each_item([this](auto & field) {
			field.template dispatch_update<System>(system_context);
		});

		controllers.for_each_item([this](auto & controller) {
			controller.template dispatch_update<System>(system_context);
		});
	}
}




