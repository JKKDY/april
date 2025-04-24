#include <vector>
#include <cassert>
#include <cstddef>

#include "debug.h"

namespace april::utils {

   
    class IndexSet {
        size_t              N;       // universe size = maxId+1
        std::vector<size_t> sparse;  // maps id to index in dense
        std::vector<size_t> dense;   // densly packed inserted IDs
    
    public:
        explicit IndexSet(size_t maxId)
          : N(maxId + 1),
            sparse(N, static_cast<size_t>(-1))
        {
            dense.reserve(N);
        }
    
        // insert (must be unique)
        void insert(size_t id) {
            // AP_ASSERT(id < N && sparse[id] == static_cast<size_t>(-1), "inserting duplicate or out‑of‑range ID");
            sparse[id] = dense.size(); // insert index into sparse
            dense.push_back(id); // insert id into dense
        }
    
        // erase then mark slot free
        void erase(size_t id) {
            // AP_ASSERT(id < N && sparse[id] < dense.size() && dense[sparse[id]] == id, "erasing non‑existent ID");
            size_t last = dense.back();
            size_t pos  = sparse[id];

            // swap last and id
            dense[pos] = last;
            sparse[last] = pos;

            dense.pop_back();
            sparse[id] = static_cast<size_t>(-1);
        }
    
        // membership test
        bool contains(size_t id) const noexcept {
            return id < N
                && sparse[id] < dense.size()
                && dense[sparse[id]] == id;
        }
    
        // iteration over live IDs
        auto begin() const noexcept { return dense.begin(); }
        auto end()   const noexcept { return dense.end();   }
        size_t size() const noexcept { return dense.size(); }
    };
}