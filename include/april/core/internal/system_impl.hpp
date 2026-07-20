#pragma once

#include <vector>
#include "april/boundaries/boundary.hpp"
#include "april/exec/policy.hpp"
#include "april/interactions/force.hpp"
#include "april/particle/attributes.hpp"
#include "../../exec/threading/scheduling.hpp"

namespace april {

	//--------------
	// EXECUTE BATCH
	//--------------
	// ToDO extract this method to a free convenience function that maps a kernel & policies to a batch for_each_pair call
	template <class SystemConfig>
	template <VectorPolicy V, container::batching::IsBatch Batch, exec::IsKernel Kernel>
	void System<SystemConfig>::execute_batch_kernel(const Batch& batch, Kernel&& kernel)  {
		using namespace april::exec;
		using namespace april::exec::internal;

		// check kernel compatibility with batch (in regard to scalar/vector mode)
		constexpr ExecutionTrait batch_traits = std::remove_cvref_t<Batch>::vector_trait;
		constexpr ExecutionMode kernel_modes = std::remove_cvref_t<Kernel>::Mode;

		constexpr ExecutionMode required_modes = required_execution_modes<batch_traits>();
		constexpr ExecutionMode valid_modes = allowed_execution_modes<V, kernel_modes>();

		static_assert((valid_modes & required_modes) == required_modes, // required modes must be a subset of valid modes
			"[APRIL] Compatibility Failure: No valid execution path found between Batch and Kernel capability sets.");

		constexpr ExecutionMode exec_mode = resolve_execution_mode<valid_modes, required_modes>();

		// run kernel on batches
		batch.template for_each_pair<exec_mode>(kernel); // TODO deduce parallel policy from batch capability and Policy P
	}



	//--------------
	// UPDATE FORCES
	//--------------
	template <class SystemConfig>
	void System<SystemConfig>::update_forces() {

		// handle pair wise (type-type) interactions
		auto update_forces_batch = [&]<container::batching::IsBatch Batch, container::batching::IsBCP BCP>(const Batch& batch, BCP && apply_bcp) {

			auto apply_batch_update =  [&] <interactions::IsForce ForceT> (const ForceT & force) APRIL_FORCE_INLINE {
				constexpr ParticleField M = ForceT::fields | ParticleField::position;

				auto kernel = [&]<bool is_packed>(auto && p1, auto && p2) APRIL_FORCE_INLINE {
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

						if constexpr (ForceT::symmetry == interactions::ForceSymmetry::Nonsymmetric) {
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
							if constexpr (ForceT::symmetry == interactions::ForceSymmetry::Antisymmetric) {
								p1.force += f_masked;
								p2.force -= f_masked;
							} else if constexpr (ForceT::symmetry == interactions::ForceSymmetry::Symmetric) {
								p1.force += f_masked;
								p2.force += f_masked;
							}
						}
					} else {
						if (r.norm_squared() > force.cutoff2()) {
							return;
						}
						if constexpr (ForceT::symmetry == interactions::ForceSymmetry::Antisymmetric) {
							vec3 f = force(p1.to_view(), p2.to_view(), r);
							p1.force += f;
							p2.force -= f;
						}  else if constexpr (ForceT::symmetry == interactions::ForceSymmetry::Symmetric) {
							vec3 f = force(p1.to_view(), p2.to_view(), r);
							p1.force += f;
							p2.force += f;
						}  else if constexpr (ForceT::symmetry == interactions::ForceSymmetry::Nonsymmetric) {
							p1.force += force(p1.to_view(), p2.to_view(), r);
							p2.force += force(p2.to_view(), p1.to_view(), -r);
						}
					}
				};

				constexpr bool force_scalar =
					(static_cast<bool>(M & ParticleField::attributes) && !particle::IsVectorizable<ParticleAttributes>)
					|| ForceT::vector_mode == exec::ExecutionMode::Scalar
					|| vector_policy == VectorPolicy::Scalar;
				constexpr VectorPolicy vp = force_scalar ? VectorPolicy::Scalar : VectorPolicy::Auto;
				execute_batch_kernel<vp>(batch, april::universal_kernel<M, ParticleField::force>(kernel));
			};

			auto [t1, t2] = batch.types;
			force_table.dispatch(t1, t2, apply_batch_update);
		};

		// handle id-id interactions
		auto update_forces_topology_batch = [&](const container::batching::IsTopologyBatch auto & batch) {
			auto apply_batch_update = [&] <interactions::IsForce ForceT> (const ForceT & force) {
				constexpr auto Read = ForceT::fields | ParticleField::position;
				constexpr auto Write = ParticleField::force;

				auto kernel = [&](auto && p1, auto && p2) APRIL_FORCE_INLINE {
					vec3 r = p2.position - p1.position;

					if constexpr (ForceT::symmetry == interactions::ForceSymmetry::Antisymmetric) {
						vec3 f = force(p1.to_view(), p2.to_view(), r);
						p1.force += f;
						p2.force -= f;
					}  else if constexpr (ForceT::symmetry == interactions::ForceSymmetry::Symmetric) {
						vec3 f = force(p1.to_view(), p2.to_view(), r);
						p1.force += f;
						p2.force += f;
					}  else if constexpr (ForceT::symmetry == interactions::ForceSymmetry::Nonsymmetric) {
						p1.force += force(p1.to_view(), p2.to_view(), r);
						p2.force += force(p2.to_view(), p1.to_view(), -r);
					}
				};
				batch.for_each_pair(april::scalar_kernel<Read, Write>(kernel));
			};

			force_table.dispatch_id(batch.representatives.first, batch.representatives.second, apply_batch_update);
		};

		for_each_particle<parallel_policy>(
			april::universal_kernel<ParticleField::force, ParticleField::force>(
				[](auto && p) { p.force = {}; } // reset forces
			)
		);

		particle_container.template invoke_for_each_interaction_batch<parallel_policy>(update_forces_batch);
		particle_container.template invoke_for_each_topology_batch<parallel_policy>(update_forces_topology_batch);
	}



	//-----------------
	// APPLY BOUNDARIES
	//-----------------
	template <class SystemConfig>
	void System<SystemConfig>::apply_boundary_conditions() {
	    particles_to_update_buffer.clear();
	    const core::Box domain_box = this->box();

		// loop through faces sequentially
	    for (DomainFace face : all_faces) {
	        const auto& compiled_boundary = boundary_table[face];

	    	// fetch work
	        std::vector<size_t> particle_ids = query_region(compiled_boundary.boundary_region);
	        if (particle_ids.empty()) continue;

	    	// partition it for parallelization
	        exec::BlockConfig config(thread_executor.num_threads());
	        auto blocks = exec::make_linear_schedule(math::Range{0, particle_ids.size()}, config);

	        // clear local buffers for this face pass
	        for (size_t i = 0; i < thread_update_buffers.size(); ++i) {
	            thread_update_buffers[i].buffer.clear();
	        }

	    	// if the particle is inside the domain
	        auto boundary_condition_inside = [&]<typename B>(const B & bc) {
	            thread_executor.execute(blocks.size(), [&](const size_t b_idx) APRIL_FORCE_INLINE {
	            	APRIL_ASSERT(exec::thread_index() < thread_executor.num_threads(),
	            		"[APRIL] exec::thread_index() must be in range [0, executor.num_threads()]. "
						"Verify that the executor sets ScopedThreadContext correctly");
	                const auto& block = blocks[b_idx];
	                auto& local_buffer = thread_update_buffers[exec::thread_index()].buffer;

	                for (size_t i = block.start; i < block.stop; ++i) {
	                    const size_t p_idx = particle_ids[i];
	                    auto p = at<B::fields>(p_idx);

	                    bc.apply(p, domain_box, face);

	                    if (compiled_boundary.topology.may_change_particle_position) {
	                        local_buffer.push_back(p_idx);
	                    }
	                }
	            });
	        };

	    	// if the particle is outside the domain
	        auto boundary_condition_outside = [&]<typename B>(const B & bc) {
	            static constexpr ParticleField detect_mask = ParticleField::position | ParticleField::old_position;

	            thread_executor.execute(blocks.size(), [&](const size_t b_idx) APRIL_FORCE_INLINE {
	                const auto& block = blocks[b_idx];
	                auto& local_buffer = thread_update_buffers[exec::thread_index()].buffer;

	                for (size_t i = block.start; i < block.stop; ++i) {
	                    const size_t p_idx = particle_ids[i];
	                    auto particle = at<B::fields | detect_mask>(p_idx);

	                    const int ax = boundary::axis_of_face(face);
	                    const vec3 diff = particle.position - particle.old_position;

	                    if (std::abs(diff[ax]) < 1e-12) continue;

	                	// we check if the particle has crossed the current face in the last time step
	                	// we solve for the intersection point via solving for t in y = t * diff + p where
	                	// diff is the displacement, p is the particles starting position and y is the face coordinate
	                    const double y = diff[ax] < 0 ? domain_box.min[ax] : domain_box.max[ax];
	                    const double t = (y - particle.old_position[ax]) / diff[ax];
	                    const vec3 intersection = t * diff + particle.old_position;

	                	// check if the intersection is part of the face
	                    auto [ax1, ax2] = boundary::non_face_axis(face);
	                    if (domain_box.max[ax1] >= intersection[ax1] && domain_box.min[ax1] <= intersection[ax1] &&
	                        domain_box.max[ax2] >= intersection[ax2] && domain_box.min[ax2] <= intersection[ax2]) {

	                        bc.apply(particle, domain_box, face);

	                        if (compiled_boundary.topology.may_change_particle_position) {
	                            local_buffer.push_back(p_idx);
	                        }
	                    }
	                }
	            });
	        };

	        if (compiled_boundary.topology.boundary_thickness >= 0) { // >0 implies the boundary region only applies to inside
	            compiled_boundary.dispatch(boundary_condition_inside);
	        } else {
	            compiled_boundary.dispatch(boundary_condition_outside);
	        }

	        // merge local buffers sequentially
	    	if (compiled_boundary.topology.may_change_particle_position) {
	    		for (size_t i = 0; i < thread_update_buffers.size(); ++i) {
	    			auto& local_buf = thread_update_buffers[i].buffer;
	    			if (local_buf.empty()) continue;

	    			particles_to_update_buffer.insert(
						particles_to_update_buffer.end(),
						local_buf.begin(),
						local_buf.end()
					);
	    		}
	    	}
	    }

		// invoke structure update on the container
	    if (!particles_to_update_buffer.empty()) {
	        particle_container.invoke_notify_moved(particles_to_update_buffer);
	    }
	}


	//------------------
	// APPLY CONTROLLERS
	//------------------
	template <class SystemConfig>
	void System<SystemConfig>::apply_controllers() {
		controllers.for_each_item([this](auto & controller) {
			if (controller.should_trigger(trig_context)) {
				controller.dispatch_apply(system_context);
			}
		});
	}


	//-------------
	// APPLY FIELDS
	//-------------
	template <class SystemConfig>
	void System<SystemConfig>::apply_force_fields() {
		fields.for_each_item([&]<typename F>(F & field) {
			for_each_particle<parallel_policy>(scalar_kernel<F::fields, ParticleField::force>(
				[&](auto && p) {
					field.dispatch_apply(p);
				})
			);
		});
	}


	//-------
	// UPDATE
	//-------
	template <class SystemConfig>
	void System<SystemConfig>::update_all_components() {
		fields.for_each_item([this](auto & field) {
			field.template dispatch_update<System>(system_context);
		});

		controllers.for_each_item([this](auto & controller) {
			controller.template dispatch_update<System>(system_context);
		});
	}
}









