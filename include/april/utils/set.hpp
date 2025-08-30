#pragma once
#include <vector>
#include <cstddef>

#include "april/utils/debug.h"

namespace april::utils {

   
    template <std::unsigned_integral T> class IndexSet {
    public:
        size_t N; // universe size = maxId+1
        std::vector<T> sparse; // maps id to index in dense
        std::vector<T> dense;  // densely packed inserted IDs

        explicit IndexSet(const T maxId = 0):
            N(maxId + 1), sparse(N, static_cast<T>(-1))
        {
            dense.reserve(N);
        }

        void set_capacity(const T maxId) {
            N = maxId + 1;
            sparse = std::vector<T>(N, static_cast<T>(-1));
            dense.reserve(N);
        }
    
        void insert(const T id) {
            AP_ASSERT(id < N && sparse[id] == static_cast<T>(-1), "inserting duplicate or out‑of‑range ID");
            sparse[id] = static_cast<T>(dense.size()); // insert index into sparse
            dense.push_back(id); // insert id into dense
        }
    
        // erase then mark slot free
        void erase(const T id) {
            AP_ASSERT(id < N && sparse[id] < dense.size() && dense[sparse[id]] == id, "erasing non‑existent ID");
            const T last = dense.back();
            const T pos  = sparse[id];

            // swap last and id
            dense[pos] = last;
            sparse[last] = pos;

            dense.pop_back();
            sparse[id] = static_cast<T>(-1);
        }
    
        // membership test
        [[nodiscard]] bool contains(const T id) const noexcept {
            return id < N
                && sparse[id] < dense.size()
                && dense[sparse[id]] == id;
        }

        // iteration over live IDs
        [[nodiscard]] auto begin() const noexcept { return dense.begin(); }
        [[nodiscard]] auto end()   const noexcept { return dense.end(); }

        [[nodiscard]] T operator[] (const size_t i) const noexcept {return dense[i]; }
        [[nodiscard]] size_t size() const noexcept { return dense.size(); }
    };
}