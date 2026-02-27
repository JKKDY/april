#pragma once
#include <utility>

#include "april/exec/policy.hpp"
#include "april/particle/packed_access.hpp"

namespace april {
    namespace exec::internal {
        // check that we can call Func with the provided args with or without a template bool
        template<typename Func, typename... Args>
        concept IsKernelInvocable =
        requires(const Func& f, Args&&... args) {
            f.template operator()<(particle::IsPackedParticleAccessor<std::remove_cvref_t<Args>> || ...)>(std::forward<Args>(args)...);
        } ||
        std::invocable<const Func&, Args...>;


        template<ExecutionMode M, typename Func>
        struct KernelWrapper {
            const Func func;
            static constexpr auto Mode = M;

            template<typename... Args>
            requires IsKernelInvocable<Func, Args...>
            decltype(auto) operator()(Args&&... args) const {
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
        template<typename F> using ScalarKernel    = KernelWrapper<ExecutionMode::Scalar, std::remove_cvref_t<F>>;
        template<typename F> using VectorKernel    = KernelWrapper<ExecutionMode::Vector, std::remove_cvref_t<F>>;
        template<typename F> using UniversalKernel = KernelWrapper<ExecutionMode::Hybrid, std::remove_cvref_t<F>>;


        // Trait to identify if a type is one of our specific wrappers
        template<typename T> struct is_kernel_wrapper : std::false_type {};

        template<ExecutionMode M, typename F>
        struct is_kernel_wrapper<KernelWrapper<M, F>> : std::true_type {};
    } // namespace exec::internal


    namespace exec {
        // Kernel concept
        template<typename T>
        concept IsKernel = internal::is_kernel_wrapper<std::remove_cvref_t<T>>::value;
    }



    template<typename F> auto scalar_kernel(F&& f) {
        return exec::internal::ScalarKernel<F>{std::forward<F>(f)};
    }
    template<typename F> auto vector_kernel(F&& f) {
        return exec::internal::VectorKernel<F>{std::forward<F>(f)};
    }
    template<typename F> auto universal_kernel(F&& f) {
        return exec::internal::UniversalKernel<F>{std::forward<F>(f)};
    }
} //namespace april









