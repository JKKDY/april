#pragma once
#include <cassert>
#include <functional>
#include <utility>
#include <unordered_set>

#include "april/env/particle.h"
#include "april/forces/force.h"
#include "april/forces/no_force.h"
#include "april/boundaries/boundary.h"
#include "april/common.h"
#include "april/env/environment.h"


namespace april::force::internal {

    // internal placeholder only
    struct ForceSentinel : Force {
        ForceSentinel() : Force(-1.0) {}

        vec3 operator()(const env::internal::Particle&, const env::internal::Particle&, const vec3&) const noexcept {
            assert(false && "NullForce should never be executed");
        }
        [[nodiscard]] ForceSentinel mix(ForceSentinel const&) const { return {}; }
    };


    template<class... Fs>
    struct VariantType {
        // 1. Disallow the internal sentinel type in user packs
        static_assert((!std::is_same_v<ForceSentinel, Fs> && ...),
                      "ForceSentinel must NOT appear in ForcePack (internal sentinel only).");

        // 2. Detect whether NoForce is already supplied
        static constexpr bool has_no_force = (std::is_same_v<NoForce, Fs> || ...);

        // 3. Compute the variant type
        using type = std::conditional_t<
            has_no_force,
            std::variant<ForceSentinel, Fs...>,           // user already included NoForce
            std::variant<ForceSentinel, Fs..., NoForce>   // append it
        >;
    };

    // Convenience alias
    template<class... Fs>
    using VariantType_t = typename VariantType<Fs...>::type;



    // TODO template just on IsForce... Fs instead of the entire environment!
    template<class Env> class ForceTable;

    template<IsForce... Fs, boundary::IsBoundary... BCs, controller::IsController...Cs, field::IsField...FFs>
    class ForceTable<env::Environment<ForcePack<Fs...>, boundary::BoundaryPack<BCs...>, controller::ControllerPack<Cs...>, field::FieldPack<FFs...>>> {
        // big variant used internally by the manager. Contains NullForce (internal placeholder) & NoForce
        using internal_force_variant_t = VariantType_t<Fs...>;

        // small variant specified by the user
        using user_force_variant_t = std::variant<Fs...>;
    public:

		ForceTable() = default;

        void build(std::vector<InteractionInfo<user_force_variant_t>>& interaction_infos,
                   const std::unordered_map<env::ParticleType, env::internal::ParticleType>& usr_types_to_impl_types,
                   const std::unordered_map<env::ParticleID, env::internal::ParticleID>& usr_ids_to_impl_ids
        ) {
            auto [type_infos, id_infos] = partition_infos(interaction_infos);
            build_type_forces(type_infos, usr_types_to_impl_types);
            build_id_forces(id_infos, usr_ids_to_impl_ids);
            validate_force_tables();
            compute_max_cutoff();
        }


        [[nodiscard]] vec3 evaluate(const env::internal::Particle& p1, const env::internal::Particle& p2) const {
            return evaluate(p1, p2, p2.position - p1.position); // dist vector points from p1 to p2
        }


        [[nodiscard]] vec3 evaluate(const env::internal::Particle& p1, const env::internal::Particle& p2, const vec3& r) const {
            auto & tF = get_type_force(p1.type, p2.type);
            vec3 force = std::visit([&](auto const& f){ return f(p1,p2,r); }, tF);

            // check if both particles even have any individual interactions defined for them
            if (p1.id < n_ids && p2.id < n_ids) {
                auto & iF = get_id_force(p1.id, p2.id);
                // TODO replace std::visit with a thunk
                force += std::visit([&](auto const& f){ return f(p1,p2,r); }, iF);
            }

            return force;
        }

        [[nodiscard]] double get_max_cutoff() const {
            return max_cutoff;
        }


    private:
        std::vector<internal_force_variant_t> inter_type_forces; // Forces between different particle types (e.g. type A <-> type B)
        std::vector<internal_force_variant_t> intra_particle_forces; // Forces between specific particle instances (by ID e.g. id1 <-> id2)

        size_t n_types{};
        size_t n_ids{};

        double max_cutoff = 0;

        [[nodiscard]] size_t type_index(const size_t a, const size_t b) const noexcept{
            return n_types * a + b;
        }

        [[nodiscard]] size_t id_index(const size_t a, const size_t b) const noexcept{
            return n_ids * a + b;
        }

        internal_force_variant_t & get_type_force(const size_t a, const size_t b) noexcept {
            return inter_type_forces[type_index(a, b)];
        }

        internal_force_variant_t & get_id_force(const size_t a, const size_t b) noexcept {
            return intra_particle_forces[id_index(a, b)];
        }

        const internal_force_variant_t& get_type_force(const size_t a, const size_t b) const noexcept {
            return inter_type_forces[type_index(a,b)];
        }

        const internal_force_variant_t& get_id_force(const size_t a, const size_t b) const noexcept {
            return intra_particle_forces[id_index(a,b)];
        }

        auto partition_infos(std::vector<InteractionInfo<user_force_variant_t>>& infos) {
            // partition type vs id
            const auto it = std::partition(infos.begin(), infos.end(),
                [](const auto& info){ return info.pair_contains_types; });

            // contains all interaction infos for particle types
            std::vector<InteractionInfo<user_force_variant_t>> type_infos{
                std::make_move_iterator(infos.begin()),
                std::make_move_iterator(it)
            };

            // contains all interaction infos for particle (id) pairs
            std::vector<InteractionInfo<user_force_variant_t>> id_infos{
                std::make_move_iterator(it),
                std::make_move_iterator(infos.end())
            };
            return std::pair{std::move(type_infos), std::move(id_infos)};
        }

        void build_type_forces(std::vector<InteractionInfo<user_force_variant_t>>& type_infos,
                       const std::unordered_map<env::ParticleType, env::internal::ParticleType>& type_map)
        {
            // helper lambda (local) to widen info-variant -> manager-variant
            auto widen = []<typename SmallV>(SmallV&& small) -> internal_force_variant_t {
                return internal_force_variant_t{std::forward<SmallV>(small)};
            };

            // collect unique particle types to define types map size (implementation types are dense [0, N-1])
            std::unordered_set<env::ParticleType> particle_types;
            for (auto& x : type_infos) {
                particle_types.insert(type_map.at(x.key_pair.first));
                particle_types.insert(type_map.at(x.key_pair.second));
            }

            n_types = particle_types.size();
            inter_type_forces.resize(n_types * n_types);

            // insert type forces into map & apply usr mappings
            for (auto& x : type_infos) {
                const auto a = type_map.at(x.key_pair.first);
                const auto b = type_map.at(x.key_pair.second);
                inter_type_forces[type_index(a, b)] = std::visit(widen, x.force);
                inter_type_forces[type_index(b, a)] = inter_type_forces[type_index(a, b)];
            }

            mix_missing_type_pairs();
        }

        void mix_missing_type_pairs() {
            //  mix missing type pairs from diagonals
            for (size_t a = 0; a < n_types; ++a) {
                for (size_t b = 0; b < n_types; ++b) {
                    auto& f = get_type_force(a, b);
                    if (a == b || !std::holds_alternative<ForceSentinel>(f)) continue;

                    auto& fa = get_type_force(a, a);
                    auto& fb = get_type_force(b, b);

                    auto force = std::visit([]<typename F1, typename F2>(F1 const& A, F2 const& B) -> internal_force_variant_t {
                            if constexpr (std::same_as<F1, F2>)
                                return A.mix(B);
                            else
                                throw std::invalid_argument("Cannot mix different force types");
                        },
                        fa, fb);

                    inter_type_forces[type_index(a, b)] = force;
                    inter_type_forces[type_index(b, a)] = force;
                }
            }
        }

        void build_id_forces(std::vector<InteractionInfo<user_force_variant_t>>& id_infos,
                     const std::unordered_map<env::ParticleID, env::internal::ParticleID>& id_map)
        {
            // helper lambda (local) to widen info-variant -> manager-variant
            auto widen = []<typename SmallV>(SmallV&& small) -> internal_force_variant_t {
                return internal_force_variant_t{std::forward<SmallV>(small)};
            };

            // collect particle ids to define ids map size (implementation ids are dense [0, M-1])
            std::unordered_set<env::ParticleID> ids;
            for (auto& x : id_infos) {
                ids.insert(id_map.at(x.key_pair.first));
                ids.insert(id_map.at(x.key_pair.second));
            }

            n_ids = ids.size();
            intra_particle_forces.resize(n_ids * n_ids);

            // insert id forces into map & apply usr mappings
            for (auto& x : id_infos) {
                const auto a = id_map.at(x.key_pair.first);
                const auto b = id_map.at(x.key_pair.second);
                intra_particle_forces[id_index(a, b)] = std::visit(widen, x.force);
                intra_particle_forces[id_index(b, a)] = intra_particle_forces[id_index(a, b)];
            }

            fill_missing_id_interactions();
        }

        void fill_missing_id_interactions() {
            // Fill undefined id interactions with no forces
            for (size_t a = 0; a < n_ids; a++) {
                for (size_t b = 0; b < n_ids; b++) {
                    auto & v = intra_particle_forces[id_index(a, b)];
                    if (a != b && std::holds_alternative<ForceSentinel>(v)) {
                        v = NoForce();
                    }
                }
            }
        }

        void validate_force_tables() const {
            #ifndef NDEBUG
            for (size_t i = 0; i < n_types; ++i)
                for (size_t j = 0; j < n_types; ++j)
                    AP_ASSERT(!std::holds_alternative<ForceSentinel>(inter_type_forces[type_index(i, j)]),
                              "inter_type_forces should not contain ForceSentinel");

            for (size_t i = 0; i < n_ids; ++i)
                for (size_t j = 0; j < n_ids; ++j) {
                    auto& v = intra_particle_forces[id_index(i, j)];
                    if (i == j)
                        AP_ASSERT(std::holds_alternative<ForceSentinel>(v),
                                  "intra_particle_forces should contain ForceSentinel for identical ids");
                    else
                        AP_ASSERT(!std::holds_alternative<ForceSentinel>(v),
                                  "intra_particle_forces should not contain ForceSentinel for differing ids");
                }
            #endif
        }


        void compute_max_cutoff() noexcept {
            // get the max cutoff distance
            auto cutoff_of = [](auto const& v){
                return std::visit([](auto const& f){ return f.cutoff; }, v);
            };

            max_cutoff = 0.0;

            // 1. Scan type-based forces
            for (auto const& v : inter_type_forces)
                max_cutoff = std::max(max_cutoff, cutoff_of(v));

            // 2. Scan id-based forces (ignore sentinels)
            for (auto const& v : intra_particle_forces)
                if (!std::holds_alternative<ForceSentinel>(v))
                    max_cutoff = std::max(max_cutoff, cutoff_of(v));
        }
    };
} // namespace april::env::impl



