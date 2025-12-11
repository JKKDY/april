#pragma once
#include <functional>
#include <set>
#include <unordered_set>

#include "april/forces/force.h"
#include "april/forces/noforce.h"


namespace april::force::internal {

    struct InteractionProp {
        double cutoff = 0.0;
        bool is_active = false;
        std::vector<std::pair<env::ParticleType, env::ParticleType>> used_by_types;
        std::vector<std::pair<env::ParticleID, env::ParticleID>> used_by_ids;
    };

    struct InteractionSchema {
        const std::vector<env::ParticleType> types;
        const std::vector<env::ParticleID> ids;

        const std::vector<InteractionProp> interactions;

        const std::vector<size_t> type_interaction_matrix; // i * types.size() + j -> index into interactions
        const std::vector<size_t> id_interaction_matrix; // i * id.size() + j -> index into interactions
    };

    template<IsForceVariant ForceVariant>
    class ForceTable {
        using Type_Interaction = TypeInteraction<ForceVariant>;
        using Id_Interaction = IdInteraction<ForceVariant>;
        using IdMap = std::unordered_map<env::ParticleID, env::ParticleID>;
        using TypeMap = std::unordered_map<env::ParticleType, env::ParticleType>;
    public:

        ForceTable(
            std::vector<Type_Interaction> type_interactions,
            std::vector<Id_Interaction> id_interactions,
            const TypeMap & usr_types_to_impl_types,
            const IdMap & usr_ids_to_impl_ids
        ) {
            build_type_forces(type_interactions, usr_types_to_impl_types);
            build_id_forces(id_interactions, usr_ids_to_impl_ids);
            validate_force_tables();
        }

        InteractionSchema generate_schema() {

            // std::vector<env::ParticleType> types(n_types);
            // std::vector<env::ParticleID> ids(n_ids);
            //
            // for (int i = 0; i < n_types; ++i) types[i] = i;
            // for (int i = 0; i < n_ids; ++i) ids[i] = i;
            //
            // InteractionProp get_properties = [](auto const& v){
            //     return std::visit([](auto const& f) {
            //         InteractionProp prop;
            //         prop.cutoff = f.cutoff();
            //     }, v);
            // };
            //
            //
            //
            //
            // return InteractionSchema{
            //     .types = types,
            //     .ids = ids
            // };
        }


        template<typename Func>
        void dispatch(const env::ParticleType t1, const env::ParticleType t2, Func && func) const {
            const auto & variant = get_type_force(t1, t2);
            std::visit([&]<IsForce F>(const F & f) -> void {
                if constexpr (!std::same_as<F, ForceSentinel> && !std::same_as<F, NoForce>) {
                    func(f);
                }
            }, variant);
        }

        template<typename Func>
        void dispatch_id(const env::ParticleID id1, const env::ParticleID id2, Func && func) const {
            const auto & variant = get_id_force(id1, id2);
            std::visit([&]<IsForce F>(const F & f) -> void {
                if constexpr (!std::same_as<F, ForceSentinel> && !std::same_as<F, NoForce>) {
                    func(f);
                }
            }, variant);
        }


        [[nodiscard]] bool has_id_force(const env::ParticleID a, const env::ParticleID b) const noexcept{
            return a < n_ids && b < n_ids;
        }

        ForceVariant & get_type_force(const env::ParticleType a, const env::ParticleType b) noexcept {
            return type_forces[type_index(a, b)];
        }

        ForceVariant & get_id_force(const env::ParticleID a, const env::ParticleID b) noexcept {
            return id_forces[id_index(a, b)];
        }

        const ForceVariant& get_type_force(const env::ParticleType a, const env::ParticleType b) const noexcept {
            return type_forces[type_index(a,b)];
        }

        const ForceVariant& get_id_force(const env::ParticleID a, const env::ParticleID b) const noexcept {
            return id_forces[id_index(a,b)];
        }


    private:
        std::vector<ForceVariant> type_forces; // Forces between different particle types (e.g. type A <-> type B)
        std::vector<ForceVariant> id_forces; // Forces between specific particle instances (by ID e.g. id1 <-> id2)
        size_t n_types{};
        size_t n_ids{};


        double max_cutoff = 0;

        [[nodiscard]] size_t type_index(const env::ParticleType a, const env::ParticleType b) const noexcept{
            return n_types * a + b;
        }

        [[nodiscard]] size_t id_index(const env::ParticleID a, const env::ParticleID b) const noexcept{
            return n_ids * a + b;
        }



        void build_type_forces(std::vector<Type_Interaction>& type_infos, const TypeMap & type_map)
        {
            // collect unique particle types to define types map size (implementation types are dense [0, N-1])
            std::unordered_set<env::ParticleType> particle_types;
            for (auto& x : type_infos) {
                particle_types.insert(type_map.at(x.type1));
                particle_types.insert(type_map.at(x.type2));
            }

            n_types = particle_types.size();
            type_forces.resize(n_types * n_types);

            // insert type forces into map & apply usr mappings
            for (auto& x : type_infos) {
                const auto a = type_map.at(x.type1);
                const auto b = type_map.at(x.type2);
                type_forces[type_index(a, b)] = x.force;
                type_forces[type_index(b, a)] = x.force;
            }

            //  mix missing type pairs from diagonals
            for (size_t a = 0; a < n_types; ++a) {
                for (size_t b = 0; b < n_types; ++b) {
                    auto& f = get_type_force(a, b);
                    if (a == b || !std::holds_alternative<ForceSentinel>(f)) continue;

                    auto& fa = get_type_force(a, a);
                    auto& fb = get_type_force(b, b);

                    auto force = std::visit([]<typename F1, typename F2>(F1 const& A, F2 const& B) -> ForceVariant {
                            if constexpr (std::same_as<F1, F2>)
                                return A.mix(B);
                            else
                                throw std::invalid_argument("Cannot mix different force types");
                        },
                        fa, fb);

                    type_forces[type_index(a, b)] = force;
                    type_forces[type_index(b, a)] = force;
                }
            }
        }


        void build_id_forces(std::vector<Id_Interaction>& id_infos, const IdMap & id_map)
        {
            // collect particle ids to define ids map size (implementation ids are dense [0, M-1])
            std::unordered_set<env::ParticleID> ids;
            for (auto& x : id_infos) {
                ids.insert(id_map.at(x.id1));
                ids.insert(id_map.at(x.id2));
            }

            n_ids = ids.size();
            id_forces.resize(n_ids * n_ids);

            // insert id forces into map & apply usr mappings
            for (auto& x : id_infos) {
                const auto a = id_map.at(x.id1);
                const auto b = id_map.at(x.id2);
                id_forces[id_index(a, b)] = x.force;
                id_forces[id_index(b, a)] = x.force;
            }

            // Fill undefined id interactions with no forces
            for (size_t a = 0; a < n_ids; a++) {
                for (size_t b = 0; b < n_ids; b++) {
                    auto & v = id_forces[id_index(a, b)];
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
                    AP_ASSERT(!std::holds_alternative<ForceSentinel>(type_forces[type_index(i, j)]),
                              "inter_type_forces should not contain ForceSentinel");

            for (size_t i = 0; i < n_ids; ++i)
                for (size_t j = 0; j < n_ids; ++j) {
                    auto& v = id_forces[id_index(i, j)];
                    if (i == j)
                        AP_ASSERT(std::holds_alternative<ForceSentinel>(v),
                                  "intra_particle_forces should contain ForceSentinel for identical ids");
                    else
                        AP_ASSERT(!std::holds_alternative<ForceSentinel>(v),
                                  "intra_particle_forces should not contain ForceSentinel for differing ids");
            }
            #endif
        }

    };
} // namespace april::env::impl



