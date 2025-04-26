
#include <cstdint>
#include <vector>
#include <utility> 
#include <limits>
#include <stdexcept>
#include <concepts>
#include <unordered_map>
#include <unordered_set>

#include "debug.h"

namespace april::utils::impl {

    static inline uint64_t splitmix64(uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x  = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x  = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        x ^= (x >> 31);
        return x;
    }

    template<typename KeyT> struct SymmetricHash {
        size_t operator()(const std::pair<KeyT, KeyT>& p) const noexcept {
            auto [a, b] = p;
            if (a > b) std::swap(a, b);
            uint64_t packed = (uint64_t(a) << 32) | uint64_t(b);
            return splitmix64(packed);
        }
    };

    template<typename KeyT> struct SymmetricEqual {
        bool operator()(const std::pair<KeyT, KeyT>& a, const std::pair<KeyT, KeyT>& b) const noexcept {
            return (a.first == b.first && a.second == b.second) || (a.first == b.second && a.second == b.first);
        }
    };


    template<std::unsigned_integral KeyT> bool keys_are_unique(const std::vector<std::pair<KeyT, KeyT>>& keys) {
        std::unordered_set<std::pair<KeyT, KeyT>, SymmetricHash<KeyT>, SymmetricEqual<KeyT>> seen;
        seen.reserve(keys.size());

        for (const auto& p : keys)
            if (!seen.insert(p).second)
                return false;
        return true;
    }


    template<typename T, std::unsigned_integral KeyT = size_t> class UnorderedMap {
        using Ptr = std::unique_ptr<T>;
    public:

        void build(const std::vector<std::pair<KeyT, KeyT>> & keys, std::vector<Ptr>&& values) {
            if (keys.size() != values.size())
                throw std::invalid_argument("keys/values size mismatch");

            if (!keys_are_unique(keys))
                throw std::invalid_argument("keys are not unique; duplicate key pairs found");

                
            for (size_t i = 0; i < keys.size(); i++) {
                auto [a,b] = keys[i];
                map[{a,b}] = std::move(values[i]);
            }
        }

        T* get(KeyT a, KeyT b) const noexcept {
            auto it = map.find({a, b});
            if (it != map.end()) {
                return it->second.get(); // return T* from the unique_ptr
            } else {
                return nullptr;
            }
        }        

        bool empty() const noexcept {
            return map.empty();
        }
    private:
        std::unordered_map<std::pair<KeyT, KeyT>, std::unique_ptr<T>, SymmetricHash<KeyT>, SymmetricEqual<KeyT>> map;
    };


    
    template<typename T, std::unsigned_integral KeyT = size_t> class DensePairMap {
        using Ptr = std::unique_ptr<T>;
    public:
    
        void build(const std::vector<std::pair<KeyT, KeyT>> & keys, std::vector<Ptr>&& values) {
            if (keys.size() != values.size())
                throw std::invalid_argument("keys/values size mismatch");
            
            if (!keys_are_unique(keys))
                throw std::invalid_argument("keys are not unique; duplicate key pairs found");

            // move in ownership
            storage = std::move(values);
            
            // allocate map
            for (const auto& [a, b] : keys) {
                if (a + 1 > N) N = a + 1;
                if (b + 1 > N) N = b + 1;
            }
            map.assign(N * N, nullptr);

            // fill map with forces
            for (size_t i = 0; i < keys.size(); ++i) {
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