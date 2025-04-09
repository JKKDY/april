#pragma once
#include <array>

#include "april/utils/Vec3.hpp"



namespace april {
	using vec3 = utils::Vec3;
	using int3 = std::array<int, 3>;
    using uint3 = std::array<uint32_t, 3>;
    using int3 = std::array<int, 3>;



    // Function taken from boost::hash_combine.
    //See: https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine.
    inline size_t hash_combine(const size_t lhs, const size_t rhs) {
        return lhs ^ (rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2));
    }

    struct Int3Hasher {
        size_t operator()(const std::array<int, 3>& arr) const {
            std::size_t seed = 0x123456789;  // initial value
            for (const int x : arr) {
                seed = hash_combine(seed, std::hash<size_t>{}(x));
            }
            return seed;
        }
    };
} // namespace april