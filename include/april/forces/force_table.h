#pragma once
#include <functional>
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

        const std::vector<size_t> type_interaction_matrix; // i * types.size() + j -> index into interactions
        const std::vector<size_t> id_interaction_matrix; // i * id.size() + j -> index into interactions

        const std::vector<InteractionProp> interactions;

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

        [[nodiscard]] InteractionSchema generate_schema() const {
            // helper to extract properties of a given force type
            auto get_properties = [](auto const& v) -> InteractionProp {
                return std::visit([]<typename F>(F const& f) -> InteractionProp {
                    using T = std::decay_t<F>;

                    InteractionProp prop;
                    prop.cutoff = f.cutoff();

                    if constexpr (std::is_same_v<T, NoForce>) {
                        prop.is_active = false;
                    } else {
                        prop.is_active = true;
                    }

                    return prop;
                }, v);
            };

            // helper to check if two variants cary the same force (regarding type as well as parameters)
            auto is_equal = [&](const ForceVariant& a, const ForceVariant& b) {
                if (a.index() != b.index()) return false;
                return std::visit([&]<typename A>(const A& val_a) {
                    using T = std::decay_t<A>;
                    const T& val_b = std::get<T>(b);
                    return val_a.equals(val_b);
                }, a);
            };

            // gather all types and ids in ascending order (types and ids are dense in [0,...])
            std::vector<env::ParticleType> types(n_types);
            std::vector<env::ParticleID> ids(n_ids);

            for (size_t i = 0; i < n_types; ++i) types[i] = i;
            for (size_t i = 0; i < n_ids; ++i) ids[i] = i;

            // first we gather all forces (type forces and id forces)
            std::vector<ForceVariant> all_forces;
            all_forces.reserve(type_forces.size() + id_forces.size());

            for (const auto & force : type_forces) all_forces.push_back(force);
            for (const auto & force : id_forces) all_forces.push_back(force);

            // and we also create a corresponding properties vector for every force
            std::vector<InteractionProp> all_force_props;
            all_force_props.reserve(type_forces.size() + id_forces.size());
            for (const auto & force : all_forces) all_force_props.push_back(get_properties(force));

            // loop through all possible type pairs and register them in the properties of their interacting force
            for (env::ParticleType i = 0; i < static_cast<env::ParticleType>(n_types); i++) {
                for (env::ParticleType j = 0; j < static_cast<env::ParticleType>(n_types); j++) {
                    all_force_props[type_index(i, j)].used_by_types.emplace_back(i, j);
                }
            }

            // loop through all possible (relevant) id pairs and register them in the properties of their interacting force
            for (env::ParticleID i = 0; i < static_cast<env::ParticleID>(n_ids); i++) {
                for (env::ParticleID j = i+1; j < static_cast<env::ParticleID>(n_ids); j++) {
                    all_force_props[type_forces.size() + id_index(i, j)].used_by_ids.emplace_back(i, j);
                }
            }

            // now we merge the properties of all identical forces
            // first we create a vector of unique forces and track which force in all_forces maps a force in unique_forces
            std::vector<size_t> remapping(all_forces.size());
            std::vector<ForceVariant> unique_forces;
            std::vector<InteractionProp> unique_props;

            for (size_t i = 0; i < all_forces.size(); i++) {
                const auto & current_force = all_forces[i];
                auto& current_prop  = all_force_props[i];

                // check if current_force is already contained in unique forces
                bool found = false;
                size_t found_idx = 0;

                for (size_t j = 0; j < unique_forces.size(); ++j) {
                    if (is_equal(current_force, unique_forces[j])) {
                        found = true;
                        found_idx = j;
                        break;
                    }
                }

                if (found) {
                    // current force is a duplicate of unique_forces[found_idx] -> merge
                    auto & props = unique_props[found_idx];

                    props.used_by_types.insert(props.used_by_types.end(), current_prop.used_by_types.begin(), current_prop.used_by_types.end());
                    props.used_by_ids.insert(props.used_by_ids.end(), current_prop.used_by_ids.begin(), current_prop.used_by_ids.end());

                    remapping[i] = found_idx;
                } else {
                    // current force is not in unique_forces -> create new entry
                    const size_t new_idx = unique_forces.size();
                    unique_forces.push_back(current_force);
                    unique_props.push_back(std::move(all_force_props[i]));

                    remapping[i] = new_idx;
                }
            }

            std::vector<size_t> type_interaction_matrix(n_types * n_types);
            std::vector<size_t> id_interaction_matrix(n_ids * n_ids);

            for (size_t i = 0; i < n_types * n_types; ++i) {
                type_interaction_matrix[i] = remapping[i];
            }

            for (size_t i = 0; i < n_ids * n_ids; ++i) {
                id_interaction_matrix[i] = remapping[n_types * n_types + i];
            }

            return InteractionSchema{
                .types = types,
                .ids = ids,
                .type_interaction_matrix = type_interaction_matrix,
                .id_interaction_matrix = id_interaction_matrix,
                .interactions = unique_props,
            };
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

            // insert type forces into map & apply user mappings
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



