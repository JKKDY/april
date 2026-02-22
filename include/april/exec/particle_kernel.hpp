#pragma once
#include <utility>

#include "april/exec/policy.hpp"
#include "april/particle/packed_access.hpp"

namespace april {
    namespace exec::internal {
        template<typename F, typename... Args>
        concept can_call_with_bool = requires(F f, Args... args) {
            f.template operator()<false>(std::forward<Args>(args)...);
            f.template operator()<true>(std::forward<Args>(args)...);
        };

        // Unified Dispatcher with boolean template value passing
        template<bool Value, typename F, typename... Args>
        decltype(auto) dispatch_kernel(F&& f, Args&&... args) {
            if constexpr (can_call_with_bool<F, Args...>) {
                return std::forward<F>(f).template operator()<Value>(std::forward<Args>(args)...);
            } else {
                return std::forward<F>(f)(std::forward<Args>(args)...);
            }
        }

        template<typename F, typename P>
        decltype(auto) call_with_index_or_particle(F&& func, size_t i, P&& p) {
            if constexpr (std::invocable<F, size_t, P>) {
                return func(i, std::forward<P>(p));
            } else if constexpr (std::invocable<F, P>) {
                return func(std::forward<P>(p));
            } else {
                static_assert(sizeof(P) == 0, "Kernel func cannot be called with these arguments");
            }
        }

        template<typename Func, typename... Args>
        concept IsKernelInvocable =
        requires(const Func& f, Args&&... args) {
            f.template operator()<(env::IsPackedParticleAccessor<std::remove_cvref_t<Args>> || ...)>(std::forward<Args>(args)...);
        } ||
        std::invocable<const Func&, Args...>;

        template<ExecutionMode M, typename Func>
        struct KernelWrapper {
            const Func func;
            static constexpr auto Mode = M;

            template<typename... Args>
            requires IsKernelInvocable<Func, Args...>
            decltype(auto) operator()(Args&&... args) const {
                constexpr bool packed = (env::IsPackedParticleAccessor<std::remove_cvref_t<Args>> || ...);

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

        template<typename F> using ScalarKernel    = KernelWrapper< ExecutionMode::Scalar, std::remove_cvref_t<F>>;
        template<typename F> using VectorKernel    = KernelWrapper< ExecutionMode::Vector, std::remove_cvref_t<F>>;
        template<typename F> using UniversalKernel = KernelWrapper< ExecutionMode::Scalar| ExecutionMode::Vector, std::remove_cvref_t<F>>;

        // // ---------------------
        // // SCALAR KERNEL WRAPPER
        // // ---------------------
        // template<typename Func>
        // struct ScalarKernel {
        //     const Func func;
        //     static constexpr auto capability = ExecutionMode::Scalar;
        //
        //     template<typename... Args>
        //     requires IsCallable<Func, Args...>
        //     decltype(auto) operator()(Args&&... args) const {
        //         constexpr bool has_packed = (env::IsPackedParticleAccessor<Args> || ...);
        //
        //         static_assert(!has_packed,
        //             "\n\nERROR: ScalarKernel invoked with SIMD types!\n"
        //             "You wrapped this function in 'scalar_kernel', but the simulation\n"
        //             "is trying to pass it Packed/SIMD data. Use 'universal_kernel' or `vector_kernel` instead.\n");
        //
        //         // Note: in c++26 use std::format to print the args
        //         // static_assert(std::is_invocable_v<const Func&, Args...>,
        //         //     "\n\nERROR: Kernel signature does not match provided arguments.\n");
        //
        //         return internal::dispatch_kernel<false>(func, std::forward<Args>(args)...);
        //     }
        // };
        //
        // // ---------------------
        // // PACKED KERNEL WRAPPER
        // // ---------------------
        // template<typename Func>
        // struct VectorKernel {
        //     const Func func;
        //     static constexpr auto capability = ExecutionMode::Vector;
        //
        //     template<typename... Args>
        //     requires (sizeof...(Args) > 0 && sizeof...(Args) <= 2) // Basic arity check
        //     decltype(auto) operator()(Args&&... args) const {
        //         constexpr bool has_scalar = (env::IsScalarParticleAccessor<Args> || ...); // if any arg is a packed
        //
        //         static_assert(!has_scalar,
        //          "\n\nERROR: VectorKernel invoked with Scalar types!\n"
        //          "You wrapped this function in 'vector_kernel', but the simulation\n"
        //          "is trying to pass it scalar data. Use 'universal_kernel' or `scalar_kernel` instead.\n");
        //
        //         // static_assert(std::is_invocable_v<const Func&, Args...>,
        //         //     "\n\nERROR: Kernel signature does not match provided arguments.\n");
        //
        //         return internal::dispatch_kernel<true>(func, std::forward<Args>(args)...);
        //     }
        // };
        //
        // // ------------------------
        // // UNIVERSAL KERNEL WRAPPER
        // // ------------------------
        // template<typename Func>
        // struct UniversalKernel {
        //     const Func func;
        //     static constexpr auto capability = ExecutionMode::Vector | ExecutionMode::Scalar;
        //
        //     template<typename... Args>
        //     requires (sizeof...(Args) > 0 && sizeof...(Args) <= 2) // Basic arity check
        //     decltype(auto) operator()(Args&&... args) const {
        //         // static_assert(std::is_invocable_v<const Func&, Args...>,
        //         //     "\n\nERROR: Kernel signature does not match provided arguments.\n");
        //
        //         constexpr bool has_packed = (env::IsPackedParticleAccessor<Args> || ...); // if any arg is a packed
        //         return internal::dispatch_kernel<has_packed>(func, std::forward<Args>(args)...);
        //     }
        // };

        // Trait to identify if a type is one of our specific wrappers
        template<typename T> struct is_kernel_wrapper : std::false_type {};

        template<ExecutionMode M, typename F>
        struct is_kernel_wrapper<KernelWrapper<M, F>> : std::true_type {};

        // template<typename F> struct is_kernel_wrapper<ScalarKernel<F>>    : std::true_type {};
        // template<typename F> struct is_kernel_wrapper<VectorKernel<F>>    : std::true_type {};
        // template<typename F> struct is_kernel_wrapper<UniversalKernel<F>> : std::true_type {};
        // template<typename F> struct is_kernel_wrapper<Is<F>> : std::true_type {};


    }

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
}

