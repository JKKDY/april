#pragma once

#include <bit>
#include <type_traits>
#include <utility>

#include "april/base/macros.hpp"
#include "april/exec/kernel.hpp"
#include "april/particle/access/policy.hpp"


namespace april::container::internal {

	template<typename...>
	inline constexpr bool always_false_v = false;


	// Helper to evaluate and execute standard prefetches at compile time
	template<typename T>
	APRIL_FORCE_INLINE void execute_prefetch(const T& field_data) {
		if constexpr (requires { field_data.prefetch(); }) {
			// It has a prefetch method (e.g. Vec3Ptr)
			field_data.prefetch();
		} else if constexpr (std::is_pointer_v<std::decay_t<T>>) {
			// It is a raw pointer (e.g. double*)
			APRIL_PREFETCH(field_data);
		}
		// If neither, do nothing
	}


	// Helper to evaluate and execute NTA prefetches at compile time
	template<typename T>
	APRIL_FORCE_INLINE void execute_prefetch_nta(const T& field_data) {
		if constexpr (requires { field_data.prefetch_nta(); }) {
			// It has a prefetch method (e.g. Vec3Ptr)
			field_data.prefetch_nta();
		} else if constexpr (std::is_pointer_v<std::decay_t<T>>) {
			// It is a raw pointer (e.g. double*)
			APRIL_PREFETCH_NTA(field_data);
		}
	}


	// Adapt a user kernel so it may optionally receive the particle index
	template<exec::IsKernel Kernel>
	auto adapt_indexed_kernel(Kernel&& kernel) {
		return exec::make_kernel_wrapper<Kernel>(
			[kernel = std::forward<Kernel>(kernel)]<bool is_packed>(
				size_t i,
				auto&& p
			) APRIL_FORCE_INLINE -> decltype(auto) {
				if constexpr (requires { kernel(i, std::forward<decltype(p)>(p)); }) {
					return kernel(i, std::forward<decltype(p)>(p)); // user wants index
				} else if constexpr (requires { kernel(std::forward<decltype(p)>(p)); }) {
					return kernel(std::forward<decltype(p)>(p)); // user only wants particle
				} else {
					// TODO in C++26 use std::format and introspection to print the received signature
					static_assert(false,
						"[APRIL] Kernel is malformed! It must have signature (size_t, auto&& p) or (auto&& p)."
					);
				}
			}
		);
	}

} // namespace april::container::internal


namespace april::container {

	// ----------------
	// PREFETCHING
	// ----------------
	template<IsContainerBuildConfig BuildConfiguration>
	template<ParticleField Mask>
	APRIL_FORCE_INLINE void Container<BuildConfiguration>::prefetch_particle(
		this const auto& self,
		auto... args
	) {
		// Standard prefetch (Temporal Locality = 3) - Use for outer loops
		if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>)
			internal::execute_prefetch(self.template invoke_get_field_ptr<ParticleField::force>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>)
			internal::execute_prefetch(self.template invoke_get_field_ptr<ParticleField::position>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>)
			internal::execute_prefetch(self.template invoke_get_field_ptr<ParticleField::velocity>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>)
			internal::execute_prefetch(self.template invoke_get_field_ptr<ParticleField::old_position>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>)
			internal::execute_prefetch(self.template invoke_get_field_ptr<ParticleField::mass>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::state>)
			internal::execute_prefetch(self.template invoke_get_field_ptr<ParticleField::state>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::type>)
			internal::execute_prefetch(self.template invoke_get_field_ptr<ParticleField::type>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::id>)
			internal::execute_prefetch(self.template invoke_get_field_ptr<ParticleField::id>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::attributes>)
			internal::execute_prefetch(self.template invoke_get_field_ptr<ParticleField::attributes>(args...));
	}


	template<IsContainerBuildConfig BuildConfiguration>
	template<ParticleField Mask>
	APRIL_FORCE_INLINE void Container<BuildConfiguration>::prefetch_particle_nta(
		this const auto& self,
		auto... args
	) {
		// Non-Temporal Access prefetch (Locality = 0) - Use for inner streaming loops
		if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>)
			internal::execute_prefetch_nta(self.template invoke_get_field_ptr<ParticleField::force>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>)
			internal::execute_prefetch_nta(self.template invoke_get_field_ptr<ParticleField::position>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>)
			internal::execute_prefetch_nta(self.template invoke_get_field_ptr<ParticleField::velocity>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>)
			internal::execute_prefetch_nta(self.template invoke_get_field_ptr<ParticleField::old_position>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>)
			internal::execute_prefetch_nta(self.template invoke_get_field_ptr<ParticleField::mass>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::state>)
			internal::execute_prefetch_nta(self.template invoke_get_field_ptr<ParticleField::state>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::type>)
			internal::execute_prefetch_nta(self.template invoke_get_field_ptr<ParticleField::type>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::id>)
			internal::execute_prefetch_nta(self.template invoke_get_field_ptr<ParticleField::id>(args...));

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::attributes>)
			internal::execute_prefetch_nta(self.template invoke_get_field_ptr<ParticleField::attributes>(args...));
	}



	// ------------------
	// PARTICLE ITERATION
	// ------------------
	template<IsContainerBuildConfig BuildConfiguration>
	template<ParallelPolicy P, VectorPolicy V, bool is_const, MaskPolicy MP, exec::IsKernel Kernel>
	void Container<BuildConfiguration>::invoke_iterate_range(
		this auto&& self,
		Kernel&& func,
		size_t start,
		size_t end
	) {
		using K = std::remove_cvref_t<Kernel>;

		constexpr auto mode = exec::internal::allowed_execution_modes<V, K::Modes>();
		auto kernel = internal::adapt_indexed_kernel(std::forward<Kernel>(func));

		self.template iterate_range<P, mode, is_const, MP>(kernel, start, end);
	}


	template<IsContainerBuildConfig BuildConfiguration>
	template<ParallelPolicy P, VectorPolicy V, bool is_const, exec::IsKernel Kernel>
	void Container<BuildConfiguration>::invoke_iterate_state(
		this auto&& self,
		Kernel&& func,
		const ParticleState state
	) {
		using K = std::remove_cvref_t<Kernel>;

		constexpr auto mode = exec::internal::allowed_execution_modes<V, K::Modes>();

		// Try optimized implementation, otherwise fall back to iterate_range.
		// The default assumes valid particle storage for the entire iteration range.
		if constexpr (requires { self.template iterate<P, mode, is_const>(func, state); }) {
			auto kernel = internal::adapt_indexed_kernel(std::forward<Kernel>(func));
			self.template iterate<P, mode, is_const>(kernel, state);
		} else {
			auto indexed_kernel = internal::adapt_indexed_kernel(func);
			auto state_filter = [&]<typename Part>(size_t i, Part&& p) {
				using PType = std::remove_cvref_t<Part>;

				if constexpr (particle::IsPackedParticleAccessor<PType>) {

					const auto mask = (p.state & state) != 0;
					p.mask_with(mask);
					indexed_kernel(i, p);

				} else if constexpr (particle::IsScalarParticleAccessor<PType>){

					if (self.index_is_valid(i) && static_cast<int>(p.state & state)) {
						indexed_kernel(i, std::forward<Part>(p));
					}

				} else {
					static_assert(false, "received non particle accessor");
				}
			};

			auto filtered_kernel = april::particle_kernel<
				mode,
				K::Read | ParticleField::state,
				K::Write
			>(state_filter);

			self.template iterate_range<P, mode, is_const, MaskPolicy::Enabled>(filtered_kernel, 0, self.capacity());
		}
	}


	//------------------------
	// PARTICLE DATA ACCESSORS
	//------------------------

	template <IsContainerBuildConfig BuildConfiguration>
	template <ParticleField Read, ParticleField Write, AccessType Access>
	auto Container<BuildConfiguration>::access_particle(this auto&& self, const auto... args) {
		constexpr bool is_const = std::is_const_v<std::remove_reference_t<decltype(self)>>;

		static_assert(
			!(is_const && Write != ParticleField::none),
			"APRIL ERROR: Cannot request write permissions (WriteMask != none) on a const Container. "
			"Either drop the write mask or ensure the container is mutable."
		);

		auto invoke = [&]<ParticleField F>(auto&& object) {
			if constexpr (Access == AccessType::Packed && requires {
				object.template get_field_ptr_packed<F>(args...);
			}) {
				return object.template get_field_ptr_packed<F>(args...);
			} else {
				return object.template get_field_ptr<F>(args...);
			}
		};

		auto get_field = [&]<ParticleField F>() {
			if constexpr (particle::internal::has_field_v<Write, F>) {
				return invoke.template operator()<F>(self);
			} else {
				const auto& const_self = self;
				return invoke.template operator()<F>(const_self);
			}
		};

		return particle::internal::make_particle_source<Read, Write>(get_field);
	}


	template<IsContainerBuildConfig BuildConfiguration>
	template<ParticleField Read, ParticleField Write, AccessType Access>
	auto Container<BuildConfiguration>::access_particle_id(
		this auto&& self,
		const ParticleID id
	) {
		// Implementing get_field_ptr_id is optional. The fallback converts
		// ID -> index and then delegates to access_particle.
		constexpr auto Mask = Read | Write;

		// Guard against none because shifting by countr_zero(0) would be invalid.
		if constexpr (Mask == ParticleField::none) {
			return self.template access_particle<Read, Write, Access>(
				self.invoke_id_to_index(id)
			);
		}

		// Pick the first active field to test whether get_field_ptr_id exists.
		// The function cannot be checked generically because it is a template.
		[[maybe_unused]] constexpr auto TestMask =
			static_cast<ParticleField>(
				1 << std::countr_zero(+(Read | Write))
			);

		// does get_field_ptr_id<TestMask>(id) compile? -> if yes, use specialized path
		if constexpr (requires { self.template get_field_ptr_id<TestMask>(id); }) {
			// Specialized path: direct ID access
			constexpr bool is_const = std::is_const_v<std::remove_reference_t<decltype(self)>>;

			static_assert(
				!(is_const && Write != ParticleField::none),
				"APRIL ERROR: Cannot request write permissions (WriteMask != none) "
				"on a const Container. Either drop the write mask or ensure the "
				"container is mutable."
			);

			auto invoke = [&]<ParticleField F>(auto&& object) {
				APRIL_ASSERT(self.contains_id(id), "Got invalid Id: " + std::to_string(id));

				if constexpr (Access == AccessType::Packed && requires {
					object.template get_field_ptr_id_packed<F>(id);
				}) {
					return object.template get_field_ptr_id_packed<F>(id);
				} else {
					return object.template get_field_ptr_id<F>(id);
				}
			};

			auto get_field = [&]<ParticleField F>() {
				if constexpr (particle::internal::has_field_v<Write, F>) {
					return invoke.template operator()<F>(self);
				} else {
					const auto& const_self = self;
					return invoke.template operator()<F>(const_self);
				}
			};

			return particle::internal::make_particle_source<Read, Write>(get_field);

		} else {
			// Fallback path: ID -> index -> access
			return self.template access_particle<Read, Write, Access>(
				self.invoke_id_to_index(id)
			);
		}
	}

} // namespace april::container