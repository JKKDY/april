#pragma once

#include <bit>
#include <type_traits>
#include <utility>
#include <vector>

#include "april/base/macros.hpp"
#include "april/exec/kernel.hpp"


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
					// TODO print kernel name by implementing a name demangler
					static_assert(false,
						"[APRIL] Kernel is malformed! It must have signature (size_t, auto&& p) or (auto&& p)."
					);
				}
			}
		);
	}


	// Adapt a kernel to stage packed memory references through register-backed buffers
	template<exec::IsKernel Kernel>
	auto adapt_buffered_kernel(Kernel&& kernel) {
		return exec::make_kernel_wrapper<Kernel>(
			[kernel = std::forward<Kernel>(kernel)]<bool is_packed>(size_t i,auto&& p) APRIL_FORCE_INLINE {
				using P = std::remove_cvref_t<decltype(p)>;

				if constexpr (particle::IsPackedParticleRef<P>) {
					static_assert(particle::IsPackedParticleAccessor<P>);

					auto buffer = p.load_buffer();
					auto view = buffer.to_view();

					kernel(i, view);
					buffer.update_into(p);
				} else {
					static_assert(particle::IsScalarParticleAccessor<P>);
					kernel(i, std::forward<decltype(p)>(p));
				}
			}
		);
	}


	// Compose the index and packed-buffer adaptations
	template<exec::IsKernel Kernel>
	auto adapt_iterator_kernel(Kernel&& kernel) {
		auto indexed_kernel = adapt_indexed_kernel(std::forward<Kernel>(kernel));
		return adapt_buffered_kernel(std::move(indexed_kernel));
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
	template<ParallelPolicy P, VectorPolicy V, bool is_const, exec::IsKernel Kernel>
	void Container<BuildConfiguration>::invoke_iterate_range(
		this auto&& self,
		Kernel&& func,
		size_t start,
		size_t end
	) {
		using K = std::remove_cvref_t<Kernel>;

		constexpr auto mode = exec::internal::allowed_execution_modes<V, K::Modes>();
		auto kernel = internal::adapt_iterator_kernel(std::forward<Kernel>(func));

		self.template iterate_range<P, mode, is_const>(kernel, start, end);
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
		auto kernel = internal::adapt_iterator_kernel(std::forward<Kernel>(func));

		// Try optimized implementation, otherwise fall back to iterate_range.
		// The default assumes valid particle storage for the entire iteration range.
		if constexpr (requires { self.template iterate<P, mode, is_const>(kernel, state); }) {
			self.template iterate<P, mode, is_const>(kernel, state);
		} else {
			// iterate_range performs no checks. Encountering memory that cannot be
			// interpreted as valid particle data, or that cannot legally be accessed,
			// may crash. This fallback is only safe if the container guarantees that
			// all accessed memory is valid particle storage.
			auto state_filter = [&]<typename Part>(size_t i, Part&& p) {
				using PType = std::remove_cvref_t<Part>;

				if constexpr (particle::IsPackedParticleAccessor<PType>) {
					static_assert(particle::IsPackedParticleRef<PType>);

					// const auto mask = (p.state.load() & +state) != 0;
					// if (!any(mask)) return; // if no particle is in requested state, skip this execution
					auto temp_buf = p.load_buffer();
					const auto mask = (temp_buf.state & +state) != 0;

					kernel(i, p.mask_with(mask));
				} else {
					static_assert(particle::IsScalarParticleAccessor<PType>);

					if (self.index_is_valid(i) && static_cast<int>(p.state & state)) {
						kernel(i, std::forward<Part>(p));
					}
				}
			};

			auto filtered_kernel = april::particle_kernel<
				mode,
				K::Read | ParticleField::state,
				K::Write
			>(state_filter);

			self.template iterate_range<P, mode, is_const>(filtered_kernel, 0, self.capacity());
		}
	}


	//------------------------
	// PARTICLE DATA ACCESSORS
	//------------------------
	template<IsContainerBuildConfig BuildConfiguration>
	template<ParticleField F>
	auto Container<BuildConfiguration>::invoke_get_field_ptr(this auto&& self, auto... args) {
		return self.template get_field_ptr<F>(args...);
	}


	template<IsContainerBuildConfig BuildConfiguration>
	template<ParticleField F>
	auto Container<BuildConfiguration>::invoke_get_field_ptr_id(this auto&& self, const ParticleID id) {
		APRIL_ASSERT(self.contains_id(id), "Got invalid Id: " + std::to_string(id));
		return self.template get_field_ptr_id<F>(id);
	}

	template<IsContainerBuildConfig BuildConfiguration>
	template<ParticleField Read, ParticleField Write>
	auto Container<BuildConfiguration>::access_particle(this auto&& self, const auto... args) {
		constexpr bool is_const = std::is_const_v<std::remove_reference_t<decltype(self)>>;

		static_assert(
			!(is_const && Write != ParticleField::none),
			"APRIL ERROR: Cannot request write permissions (WriteMask != none) on a const Container. "
			"Either drop the write mask or ensure the container is mutable."
		);

		particle::internal::ParticleSource<Read, Write, ParticleAttributes> src;
		constexpr auto Mask = Read | Write;

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>)
			src.force = self.template invoke_get_field_ptr<ParticleField::force>(args...);

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>)
			src.position = self.template invoke_get_field_ptr<ParticleField::position>(args...);

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>)
			src.velocity = self.template invoke_get_field_ptr<ParticleField::velocity>(args...);

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>)
			src.old_position = self.template invoke_get_field_ptr<ParticleField::old_position>(args...);

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>)
			src.mass = self.template invoke_get_field_ptr<ParticleField::mass>(args...);

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::state>)
			src.state = self.template invoke_get_field_ptr<ParticleField::state>(args...);

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::type>)
			src.type = self.template invoke_get_field_ptr<ParticleField::type>(args...);

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::id>)
			src.id = self.template invoke_get_field_ptr<ParticleField::id>(args...);

		if constexpr (particle::internal::has_field_v<Mask, ParticleField::attributes>)
			src.attributes = self.template invoke_get_field_ptr<ParticleField::attributes>(args...);

		return src;
	}


	template<IsContainerBuildConfig BuildConfiguration>
	template<ParticleField Read, ParticleField Write>
	auto Container<BuildConfiguration>::access_particle_id(
		this auto&& self,
		const ParticleID id
	) {
		// Implementing get_field_ptr_id is optional. The fallback converts
		// ID -> index and then delegates to access_particle.
		constexpr auto Mask = Read | Write;

		// Guard against none because shifting by countr_zero(0) would be invalid.
		if constexpr (Mask == ParticleField::none) {
			return particle::internal::ParticleSource<Read,Write,ParticleAttributes>{};
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

			particle::internal::ParticleSource<Read, Write, ParticleAttributes> src;

			if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>)
				src.force = self.template invoke_get_field_ptr_id<ParticleField::force>(id);

			if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>)
				src.position = self.template invoke_get_field_ptr_id<ParticleField::position>(id);

			if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>)
				src.velocity = self.template invoke_get_field_ptr_id<ParticleField::velocity>(id);

			if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>)
				src.old_position = self.template invoke_get_field_ptr_id<ParticleField::old_position>(id);

			if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>)
				src.mass = self.template invoke_get_field_ptr_id<ParticleField::mass>(id);

			if constexpr (particle::internal::has_field_v<Mask, ParticleField::state>)
				src.state = self.template invoke_get_field_ptr_id<ParticleField::state>(id);

			if constexpr (particle::internal::has_field_v<Mask, ParticleField::type>)
				src.type = self.template invoke_get_field_ptr_id<ParticleField::type>(id);

			if constexpr (particle::internal::has_field_v<Mask, ParticleField::id>)
				src.id = self.template invoke_get_field_ptr_id<ParticleField::id>(id);

			if constexpr (particle::internal::has_field_v<Mask, ParticleField::attributes>)
				src.attributes = self.template invoke_get_field_ptr_id<ParticleField::attributes>(id);

			return src;
		} else {
			// Fallback path: ID -> index -> access
			return self.template access_particle<Read, Write>(
				self.invoke_id_to_index(id)
			);
		}
	}

} // namespace april::container