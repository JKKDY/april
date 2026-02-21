#pragma once
#include <utility>

#include "april/exec/policy.hpp"
#include "april/particle/scalar_access.hpp"
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
                // forward template bool
                return std::forward<F>(f).template operator()<Value>(std::forward<Args>(args)...);
            } else {
                // no bool forwarding
                return std::forward<F>(f)(std::forward<Args>(args)...);
            }
        }

        // ---------------------
        // SCALAR KERNEL WRAPPER
        // ---------------------
        template<typename Func>
        struct ScalarKernel {
            const Func func;
            static constexpr auto capability = ExecutionMode::Scalar;

            // Unary
            template<env::IsScalarParticleAccessor P>
            requires std::is_invocable_v<const Func&, P>
            decltype(auto) operator()(P&& p) const {
                return internal::dispatch_kernel<false>(func, std::forward<P>(p));
            }

            // Binary
            template<env::IsScalarParticleAccessor P1, env::IsScalarParticleAccessor P2>
            requires std::same_as<std::remove_cvref_t<P1>, std::remove_cvref_t<P2>> &&
                     std::is_invocable_v<const Func&, P1, P2>
            decltype(auto) operator()(P1&& p1, P2&& p2) const {
                return internal::dispatch_kernel<false>(func, std::forward<P1>(p1), std::forward<P2>(p2));
            }
        };

        // ---------------------
        // PACKED KERNEL WRAPPER
        // ---------------------
        template<typename Func>
        struct VectorKernel {
            const Func func;
            static constexpr auto capability = ExecutionMode::Vector;


            // Unary
            template<env::IsPackedParticleAccessor P>
            requires std::is_invocable_v<const Func&, P>
            decltype(auto) operator()(P&& p) const {
                return internal::dispatch_kernel<true>(func, std::forward<P>(p));
            }

            // Binary
            template<env::IsPackedParticleAccessor P1, env::IsPackedParticleAccessor P2>
            requires std::same_as<std::remove_cvref_t<P1>, std::remove_cvref_t<P2>> &&
                     std::is_invocable_v<const Func&, P1, P2>
            decltype(auto) operator()(P1&& p1, P2&& p2) const {
                return internal::dispatch_kernel<true>(func, std::forward<P1>(p1), std::forward<P2>(p2));
            }
        };

        // ------------------------
        // UNIVERSAL KERNEL WRAPPER
        // ------------------------
        template<typename Func>
        struct UniversalKernel {
            const Func func;
            static constexpr auto capability = ExecutionMode::Vector | ExecutionMode::Scalar;

            // Unary
            template<env::IsAnyParticleAccessor P>
            requires std::is_invocable_v<const Func&, P>
            decltype(auto) operator()(P&& p) const {
                constexpr bool is_packed = env::IsPackedParticleAccessor<P>;
                return internal::dispatch_kernel<is_packed>(func, std::forward<P>(p));
            }

            // Binary
            template<env::IsAnyParticleAccessor P1, env::IsAnyParticleAccessor P2>
            requires std::same_as<std::remove_cvref_t<P1>, std::remove_cvref_t<P2>> &&
                     std::is_invocable_v<const Func&, P1, P2>
            decltype(auto) operator()(P1&& p1, P2&& p2) const {
                constexpr bool is_packed = env::IsPackedParticleAccessor<P1> && env::IsPackedParticleAccessor<P2>;
                return internal::dispatch_kernel<is_packed>(func, std::forward<P1>(p1), std::forward<P2>(p2));
            }
        };

        // Trait to identify if a type is one of our specific wrappers
        template<typename T> struct is_kernel_wrapper : std::false_type {};

        template<typename F> struct is_kernel_wrapper<ScalarKernel<F>>    : std::true_type {};
        template<typename F> struct is_kernel_wrapper<VectorKernel<F>>    : std::true_type {};
        template<typename F> struct is_kernel_wrapper<UniversalKernel<F>> : std::true_type {};

        // Dummy type for probing arity. Can convert into anything implicitly
        struct UniversalAny {
            template<typename T> operator T&();
            template<typename T> operator T&&();
            template<typename T> operator const T&() const;
        };
    }

    namespace exec {
        // Base Kernel Check
        template<typename T>
        concept IsKernel = internal::is_kernel_wrapper<std::remove_cvref_t<T>>::value;

        // Unary Kernel Concept
        template<typename T>
        concept IsUnaryKernel = IsKernel<T> && requires(std::remove_cvref_t<T> k) {
            k.func(std::declval<internal::UniversalAny>());
        };

        // Binary Check
        template<typename T>
        concept IsBinaryKernel = IsKernel<T> && requires(std::remove_cvref_t<T> k) {
            k.func(std::declval<internal::UniversalAny>(), std::declval<internal::UniversalAny>());
        };
    }


    // Maker functions
    template<typename Func> auto scalar_kernel(Func&& f) {
        return exec::internal::ScalarKernel<std::remove_cvref_t<Func>>{std::forward<Func>(f)};
    }
    template<typename Func> auto vector_kernel(Func&& f) {
        return exec::internal::VectorKernel<std::remove_cvref_t<Func>>{std::forward<Func>(f)};
    }
    template<typename Func> auto universal_kernel(Func&& f) {
        return exec::internal::UniversalKernel<std::remove_cvref_t<Func>>{std::forward<Func>(f)};
    }
}
