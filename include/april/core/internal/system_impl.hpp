#pragma once

#include "april/boundaries/boundary.hpp"
#include "april/exec/policy.hpp"
#include "april/forces/force.hpp"
#include "april/particle/attributes.hpp"

namespace april {

	//--------------
	// EXECUTE BATCH
	//--------------
	template <class ContainerDecl, core::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<
		ContainerDecl, Traits>
	template <ParallelPolicy P, VectorPolicy V, container::batching::IsBatch Batch, exec::IsKernel Kernel>
	void System<ContainerDecl, Traits>::execute_batch_kernel(const Batch& batch, Kernel&& kernel)  {
		using namespace april::exec::internal;
		constexpr VectorTrait batch_traits = std::remove_cvref_t<Batch>::vector_trait;
		constexpr ExecutionMode kernel_modes = std::remove_cvref_t<Kernel>::Mode;

		constexpr ExecutionMode required_modes = required_execution_modes<batch_traits>();
		constexpr ExecutionMode valid_modes = valid_execution_modes<V, kernel_modes>();

		static_assert((valid_modes & required_modes) == required_modes, // required modes must be a subset of valid modes
			"[APRIL] Compatibility Failure: No valid execution path found between Batch and Kernel capability sets.");

		constexpr ExecutionMode exec_mode = resolve_execution_mode<valid_modes, required_modes>();

		if constexpr (container::batching::IsBatchAtom<Batch>) {
			batch.template for_each_pair<P, exec_mode>(kernel);
		}
		else if constexpr (container::batching::IsBatchAtomRange<Batch>) {
			for (const auto& atom : batch) {
				atom.template for_each_pair<P, exec_mode>(kernel);
			}
		}
	}



	//--------------
	// UPDATE FORCES
	//--------------
	template <class C, core::internal::IsEnvironmentTraits Traits> requires container::IsContainerDecl<C, Traits>
	void System<C, Traits>::update_forces() {

		// batch update lambda. passed into container::for_each_interaction_batch
		auto update_batch = [&]<container::batching::IsBatch Batch, container::batching::IsBCP BCP>(const Batch& batch, BCP && apply_bcp) {

			auto apply_batch_update =  [&] <force::IsForce ForceT> (const ForceT & force) {
				constexpr ParticleField M = ForceT::fields | ParticleField::position;

				auto kernel = [&]<bool is_packed>(auto && p1, auto && p2) {
					auto diff = p2.position - p1.position;

					const auto r = [&] {
						if constexpr (std::is_same_v<std::decay_t<BCP>, container::batching::NoBatchBCP>) {
							return diff;
						} else {
							return apply_bcp(diff);
						}
					}();

					if constexpr (is_packed) {
						auto outside  = r.norm_squared() > force.cutoff2();
						if (all(outside )) return;

						if constexpr (ForceT::symmetry == force::ForceSymmetry::Nonsymmetric) {
							auto f1 = force(p1, p2, r);
							auto f2 = force(p2, p1, -r);

							p1.force += pvec3 {
								select(outside , packed(0), f1.x),
								select(outside , packed(0), f1.y),
								select(outside , packed(0), f1.z)
							};

							p2.force += pvec3 {
								select(outside , packed(0), f2.x),
								select(outside , packed(0), f2.y),
								select(outside , packed(0), f2.z)
							};
						} else {
							auto f = force(p1, p2, r);
							auto f_masked = pvec3 {
								select(outside , packed(0), f.x),
								select(outside , packed(0), f.y),
								select(outside , packed(0), f.z)
							};
							if constexpr (ForceT::symmetry == force::ForceSymmetry::Antisymmetric) {
								p1.force += f_masked;
								p2.force -= f_masked;
							} else if constexpr (ForceT::symmetry == force::ForceSymmetry::Symmetric) {
								p1.force += f_masked;
								p2.force += f_masked;
							}
						}
					} else {
						if (r.norm_squared() > force.cutoff2()) {
							return;
						}
						if constexpr (ForceT::symmetry == force::ForceSymmetry::Antisymmetric) {
							vec3 f = force(p1.to_view(), p2.to_view(), r);
							p1.force += f;
							p2.force -= f;
						}  else if constexpr (ForceT::symmetry == force::ForceSymmetry::Symmetric) {
							vec3 f = force(p1.to_view(), p2.to_view(), r);
							p1.force += f;
							p2.force += f;
						}  else if constexpr (ForceT::symmetry == force::ForceSymmetry::Nonsymmetric) {
							p1.force += force(p1.to_view(), p2.to_view(), r);
							p2.force += force(p2.to_view(), p1.to_view(), -r);
						}
					}
				};


				constexpr bool force_scalar =
					(static_cast<bool>(M & ParticleField::attributes) && !particle::IsVectorizable<ParticleAttributes>)
					|| ForceT::VectorMode == exec::internal::ExecutionMode::Scalar;
				constexpr VectorPolicy vp = force_scalar ? VectorPolicy::Scalar : VectorPolicy::Auto;
				execute_batch_kernel<ParallelPolicy::Serial, vp>(batch, april::universal_kernel<M, ParticleField::force>(kernel));
			};

			auto [t1, t2] = batch.types;
			force_table.dispatch(t1, t2, apply_batch_update);
		};

		particle_container.for_each_particle(
			april::universal_kernel<ParticleField::force, ParticleField::force>(
				[](auto && p) { p.force = {}; } // reset forces
			)
		);
		particle_container.invoke_for_each_interaction_batch(update_batch);




		// handle id interactions
		auto update_global_batch = [&](const auto & batch) {

			auto apply_batch_update = [&] <force::IsForce ForceT> (const ForceT & force) {
				constexpr ParticleField Read = ForceT::fields | ParticleField::force | ParticleField::position;

				for (const auto & [id1, id2] : batch.pairs) {
					auto && p1 = at_id<Read, ForceT::fields | ParticleField::force>(id1);
					auto && p2 = at_id<Read, ForceT::fields | ParticleField::force>(id2);

					vec3 r = p2.position - p1.position;

					if constexpr (ForceT::symmetry == force::ForceSymmetry::Antisymmetric) {
						vec3 f = force(p1.to_view(), p2.to_view(), r);
						p1.force += f;
						p2.force -= f;
					}  else if constexpr (ForceT::symmetry == force::ForceSymmetry::Symmetric) {
						vec3 f = force(p1.to_view(), p2.to_view(), r);
						p1.force += f;
						p2.force += f;
					}  else if constexpr (ForceT::symmetry == force::ForceSymmetry::Nonsymmetric) {
						p1.force += force(p1.to_view(), p2.to_view(), r);
						p2.force += force(p2.to_view(), p1.to_view(), -r);
					}
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

		for (DomainFace face : all_faces) {

			const auto& compiled_boundary = boundary_table[face];

			std::vector<size_t> particle_ids = particle_container.invoke_collect_indices_in_region(compiled_boundary.boundary_region);

			auto boundary_condition_inside = [&]<typename B>(const B & bc) {
				constexpr ParticleField M = std::decay_t<B>::fields;

				for (auto p_idx : particle_ids) {
					auto p = at<M>(p_idx);
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
					auto particle = at<M>(p_idx);

					// make sure the particle exited through the current boundary face
					// solve for intersection of the particles path with the boundary face
					// with the equation y = t * diff + p where:
					// diff is the path traveled, p is the particles starting position and y is the face

					// TODO Review for numerical robustness
					const int ax = boundary::axis_of_face(face);
					const vec3 diff = particle.position - particle.old_position;
					const double y = diff[ax] < 0 ? domain_box.min[ax] : domain_box.max[ax];
					const double t = (y - particle.old_position[ax]) / diff[ax];

					const vec3 intersection = t * diff + particle.old_position;

					// and check if that point is on the domains surface
					auto [ax1, ax2] = boundary::non_face_axis(face);
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
		fields.for_each_item([&]<typename F>(F & field) {
			for_each_particle(scalar_kernel<F::fields, ParticleField::force>(
				[&](auto && p) {
					field.dispatch_apply(p);
				})
			);
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








