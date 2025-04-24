
#include <cstdint>
#include <vector>
#include <utility> 
#include <limits>
#include <stdexcept>
#include <concepts>
#include <unordered_map>

#include "debug.h"

namespace april::utils::impl {

    template<typename T, std::integral KeyT = size_t> class UnorderedMap {
        using Ptr = std::unique_ptr<T>;
    public:

        void build(const std::vector<std::pair<KeyT, KeyT>> & keys, std::vector<Ptr>&& values) {
            size_t N = values.size();
            if (keys.size() != N)
                throw std::invalid_argument("keys/values size mismatch");
                
            for (size_t i = 0; i < N; i++) {
                map[keys[i]] = values[i];
            }
        }

        T* get(KeyT a, KeyT b) const noexcept {
            const auto it = map.find({a,b});
            return it != map.end()? *it : nullptr;
        }

        bool empty() const noexcept {
            return map.empty();
        }
    private:
        std::unordered_map<std::pair<KeyT, KeyT>, std::unique_ptr<T>> map;
    };


    
    template<typename T, std::integral KeyT = size_t> class DensePairMap {
        using Ptr = std::unique_ptr<T>;
    public:
    
        void build(const std::vector<std::pair<KeyT, KeyT>> & keys, std::vector<Ptr>&& values) {
            N = values.size();
            if (keys.size() != N)
                throw std::invalid_argument("keys/values size mismatch");

            // move in ownership
            storage = std::move(values);
            map.assign(N * N, nullptr);

            for (size_t i = 0; i < N; ++i) {
                auto [a, b] = keys[i];

                if (static_cast<size_t>(a) >= N || static_cast<size_t>(b) >= N)
                    throw std::out_of_range("key index out of range");

                T* ptr = storage[i].get();
                map[static_cast<size_t>(a) * N + static_cast<size_t>(b)] = ptr;
                map[static_cast<size_t>(b) * N + static_cast<size_t>(a)] = ptr;
            }
        }
    
        T* get(KeyT a, KeyT b) const noexcept {
            size_t ai = static_cast<size_t>(a);
            size_t bi = static_cast<size_t>(b);
            AP_ASSERT(ai < N && bi < N, "key out of range");
            return map[ai*N + bi];
        }

        size_t key_size() const noexcept {
            return N;
        }
    
    private:
        size_t N = 0;               // # of unique keys
        std::vector<T*> map;        // size N*N, raw lookup table
        std::vector<Ptr> storage;   // owns all the unique_ptr<T>
    };
} // namespace april::utils