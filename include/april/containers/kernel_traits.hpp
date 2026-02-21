#pragma once
#include <type_traits>
#include <concepts>
#include "april/base/policy.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/scalar_access.hpp"
#include "april/particle/packed_access.hpp"


namespace april::container::internal {

    static constexpr auto TestMask = ParticleField::all;
    using TestUserData = env::NoUserData;
    using TestScalarRef = env::ScalarRestrictedParticleRef<TestMask, TestUserData>;
    using TestVectorBuf = env::PackedParticleBuffer<TestMask>;

    // Check for Scalar Kernel
    template<typename Kernel>
    concept IsScalarUnaryKernel = requires(Kernel k, TestScalarRef& p) {
        requires std::invocable<Kernel, TestScalarRef&>;
    };
    template<typename Kernel>
    concept IsScalarBinaryKernel = requires(Kernel k, TestScalarRef& p) {
        requires std::invocable<Kernel, TestScalarRef&, TestScalarRef&>;
    };

    // Check Vector Kernel
    template<typename Kernel>
    concept IsVectorUnaryKernel = requires(Kernel k, TestVectorBuf& v) {
        requires std::invocable<Kernel, TestVectorBuf&>;
    };
    template<typename Kernel>
    concept IsVectorBinaryKernel = requires(Kernel k, TestVectorBuf& v) {
        requires std::invocable<Kernel, TestVectorBuf&, TestVectorBuf&>;
    };

    // Combined unary & binary checks
    // TODO use this to restrict the kernels passed to the container/system
    template<typename K>
    concept IsUnaryKernel = IsScalarUnaryKernel<K> || IsVectorUnaryKernel<K>;

    template<typename K>
    concept IsBinaryKernel = IsScalarBinaryKernel<K> || IsVectorBinaryKernel<K>;


    // get vector trait
    template<typename Kernel>
    struct KernelTraitResolver {
        using VT = april::internal::ExecutionMode;
        static constexpr VT value = [] {
            constexpr bool s = IsScalarUnaryKernel<Kernel> || IsScalarBinaryKernel<Kernel>;
            constexpr bool v = IsVectorUnaryKernel<Kernel> || IsVectorBinaryKernel<Kernel>;


            if constexpr (s && v) return VT::Scalar | VT::Vector;
            if constexpr (s)      return VT::Scalar;
            if constexpr (v)      return VT::Vector;

            static_assert(s || v, "Kernel must satisfy at least one Scalar or Vector concept.");
            return VT::Scalar;
        }();
    };

    template<typename Kernel>
    constexpr april::internal::ExecutionMode KernelTrait = KernelTraitResolver<std::remove_cvref_t<Kernel>>::value;
}


namespace april {
    namespace internal {
        // ---------------------
        // SCALAR KERNEL WRAPPER
        // ---------------------
        template<typename Func>
        struct ScalarKernel {
            Func func;

            template<env::IsScalarParticleAccessor P>
            requires std::is_invocable_v<const Func&, P>
            decltype(auto) operator()(P&& p) const {
                return func(std::forward<P>(p));
            }

            template<env::IsScalarParticleAccessor P1, env::IsScalarParticleAccessor P2>
            requires std::same_as<std::remove_cvref_t<P1>, std::remove_cvref_t<P2>> &&
                     std::is_invocable_v<const Func&, P1, P2>
            decltype(auto) operator()(P1&& p1, P2&& p2) const {
                return func(std::forward<P1>(p1), std::forward<P2>(p2));
            }
        };

        // ---------------------
        // PACKED KERNEL WRAPPER
        // ---------------------
        template<typename Func>
        struct VectorKernel {
            Func func;

            template<env::IsPackedParticleAccessor P>
            requires std::is_invocable_v<const Func&, P>
            decltype(auto) operator()(P&& p) const {
                return func(std::forward<P>(p));
            }

            template<env::IsPackedParticleAccessor P1, env::IsPackedParticleAccessor P2>
            requires std::same_as<std::remove_cvref_t<P1>, std::remove_cvref_t<P2>> &&
                     std::is_invocable_v<const Func&, P1, P2>
            decltype(auto) operator()(P1&& p1, P2&& p2) const {
                return func(std::forward<P1>(p1), std::forward<P2>(p2));
            }
        };

        // ------------------------
        // UNIVERSAL KERNEL WRAPPER
        // ------------------------
        template<typename Func>
        struct UniversalKernel {
            Func func;

            template<env::IsAnyParticleAccessor P>
            requires std::is_invocable_v<const Func&, P>
            decltype(auto) operator()(P&& p) const {
                return func(std::forward<P>(p));
            }

            template<env::IsAnyParticleAccessor P1, env::IsAnyParticleAccessor P2>
            requires std::same_as<std::remove_cvref_t<P1>, std::remove_cvref_t<P2>> &&
                     std::is_invocable_v<const Func&, P1, P2>
            decltype(auto) operator()(P1&& p1, P2&& p2) const {
                return func(std::forward<P1>(p1), std::forward<P2>(p2));
            }
        };
    }


    // Maker functions
    template<typename Func> auto scalar_kernel(Func&& f) {
        return internal::ScalarKernel<std::remove_cvref_t<Func>>{std::forward<Func>(f)};
    }
    template<typename Func> auto vector_kernel(Func&& f) {
        return internal::VectorKernel<std::remove_cvref_t<Func>>{std::forward<Func>(f)};
    }
    template<typename Func> auto universal_kernel(Func&& f) {
        return internal::UniversalKernel<std::remove_cvref_t<Func>>{std::forward<Func>(f)};
    }
}