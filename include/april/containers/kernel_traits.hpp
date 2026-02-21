#pragma once
#include <type_traits>
#include <concepts>
#include "../exec/policy.hpp"
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
    template<typename K>
    concept IsUnaryKernel = IsScalarUnaryKernel<K> || IsVectorUnaryKernel<K>;

    template<typename K>
    concept IsBinaryKernel = IsScalarBinaryKernel<K> || IsVectorBinaryKernel<K>;


    // get vector trait
    template<typename Kernel>
    struct KernelTraitResolver {
        using VT = exec::internal::ExecutionMode;
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
    constexpr exec::internal::ExecutionMode KernelTrait = KernelTraitResolver<std::remove_cvref_t<Kernel>>::value;
}

