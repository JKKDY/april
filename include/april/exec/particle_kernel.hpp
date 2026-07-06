#pragma once

#include <concepts>
#include <type_traits>
#include <utility>

#include "april/exec/policy.hpp"
#include "april/particle/packed_access.hpp"
#include "april/particle/particle_types.hpp"


namespace april {
    namespace exec::internal {
        template<typename T>
        struct is_kernel_wrapper;

        // check that we can call Func with the provided args with or without a template bool
        template<typename Func, typename... Args>
        concept IsKernelInvocable =
        requires(const Func& f, Args&&... args) {
            f.template operator()<(particle::IsPackedParticleView<std::remove_cvref_t<Args>> || ...)>(std::forward<Args>(args)...);
        } ||
        std::invocable<const Func&, Args...>;


        template<ParticleField access_fields, ParticleField update_fields, ExecutionMode exec_mode, typename Func>
        struct KernelWrapper {
            const Func func;
            static constexpr auto Mode = exec_mode;
            static constexpr auto Read = access_fields;
            static constexpr auto Write = update_fields;

            template<typename... Args>
            requires IsKernelInvocable<Func, Args...>
            APRIL_FORCE_INLINE decltype(auto) operator()(Args&&... args) const {
                constexpr bool packed = (particle::IsPackedParticleAccessor<std::remove_cvref_t<Args>> || ...);

                // Mode checks
                if constexpr (Mode == ExecutionMode::Scalar)
                    static_assert(!packed, "Error: ScalarKernel invoked with SIMD data!");
                if constexpr (Mode == ExecutionMode::Vector)
                    static_assert(packed, "Error: VectorKernel invoked with Scalar data!");

                // Dispatch
                if constexpr (requires {func.template operator()<packed>(std::forward<Args>(args)...); }) {
                    return func.template operator()<packed>(std::forward<Args>(args)...);
                } else {
                    return func(std::forward<Args>(args)...);
                }
            }
        };


        // specialization aliases
        template<ParticleField Read, ParticleField Write, typename F> using ScalarKernel =
            KernelWrapper<Read, Write, ExecutionMode::Scalar, std::remove_cvref_t<F>>;
        template<ParticleField Read, ParticleField Write, typename F> using VectorKernel =
            KernelWrapper<Read, Write, ExecutionMode::Vector, std::remove_cvref_t<F>>;
        template<ParticleField Read, ParticleField Write, typename F> using UniversalKernel =
            KernelWrapper<Read, Write, ExecutionMode::Hybrid, std::remove_cvref_t<F>>;


        // Trait to identify if a type is one of our specific wrappers
        template<typename T> struct is_kernel_wrapper : std::false_type {};

        template<ParticleField A, ParticleField U, ExecutionMode M, typename F>
        struct is_kernel_wrapper<KernelWrapper<A, U, M, F>> : std::true_type {};
    } // namespace exec::internal


    namespace exec {
        // Kernel concept
        template<typename T>
        concept IsKernel = internal::is_kernel_wrapper<std::remove_cvref_t<T>>::value;

        template<IsKernel K, typename F>
        auto make_kernel_wrapper(F&& func) {
            using KernelType = std::remove_cvref_t<K>;

            return exec::internal::KernelWrapper<
                KernelType::Read,
                KernelType::Write,
                KernelType::Mode,
                std::remove_cvref_t<F>
            >{std::forward<F>(func)};
        }
    }


    template<ParticleField Read=ParticleField::none, ParticleField Write=ParticleField::none, typename F> auto scalar_kernel(F&& f) {
        return exec::internal::ScalarKernel<Read, Write, F>{std::forward<F>(f)};
    }
    template<ParticleField Read=ParticleField::none, ParticleField Write=ParticleField::none, typename F> auto vector_kernel(F&& f) {
        return exec::internal::VectorKernel<Read, Write, F>{std::forward<F>(f)};
    }
    template<ParticleField Read=ParticleField::none, ParticleField Write=ParticleField::none, typename F> auto universal_kernel(F&& f) {
        return exec::internal::UniversalKernel<Read, Write, F>{std::forward<F>(f)};
    }
} //namespace april


