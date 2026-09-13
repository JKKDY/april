#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

#include "april/interactions/interaction.hpp"
#include "april/interactions/no_interaction.hpp"


namespace april::interaction::internal {


    struct InteractionDescriptor {
        double cutoff = 0.0;
        bool is_active = false;
        size_t arity = 2; // future stub for bonded interactions. Currently unused
        std::vector<std::pair<ParticleType, ParticleType>> used_by_types;
        std::vector<std::pair<ParticleID, ParticleID>> used_by_ids;
    };


    struct InteractionMap {
        const std::vector<ParticleType> types; // list of all types
        const std::vector<ParticleID> ids;  // list of all ids
        const std::vector<InteractionDescriptor> interactions; // list of interactions

        const std::vector<size_t> type_interaction_matrix; // i * types.size() + j -> index into interactions
        const std::vector<size_t> id_interaction_matrix; // i * id.size() + j -> index into interactions
    };


    template<IsInteractionVariant InteractionVariant>
    class InteractionTable {
        using Type_Interaction = TypeInteraction<InteractionVariant>;
        using Id_Interaction = IdInteraction<InteractionVariant>;
        using IdMap = std::unordered_map<ParticleID, ParticleID>;
        using TypeMap = std::unordered_map<ParticleType, ParticleType>;
    public:

        InteractionTable(
            std::vector<Type_Interaction> type_interactions,
            std::vector<Id_Interaction> id_interactions,
            const TypeMap & usr_types_to_impl_types,
            const IdMap & usr_ids_to_impl_ids
        ) {
            build_type_interactions(type_interactions, usr_types_to_impl_types);
            build_id_interactions(id_interactions, usr_ids_to_impl_ids);
            validate_interaction_table();
        }

        [[nodiscard]] InteractionMap generate_interaction_map() const {
            // helper to extract properties of a given interaction type
            auto get_properties = [](auto const& variant) -> InteractionDescriptor {
                return std::visit([]<typename I>(I const& interaction) -> InteractionDescriptor {
                    using T = std::decay_t<I>;

                    InteractionDescriptor prop;
                    prop.cutoff = interaction.cutoff();

                    if constexpr (std::is_same_v<T, NoInteraction>) {
                        prop.is_active = false;
                    } else {
                        prop.is_active = true;
                    }

                    return prop;
                }, variant);
            };

            // helper to check if two variants cary the same interaction (regarding type as well as parameters)
            auto is_equal = [&](const InteractionVariant& a, const InteractionVariant& b) {
                if (a.index() != b.index()) return false;
                return std::visit([&]<typename A>(const A& val_a) {
                    using T = std::decay_t<A>;
                    const T& val_b = std::get<T>(b);
                    return val_a.equals(val_b);
                }, a);
            };

            // gather all types and ids in ascending order (types and ids are dense in [0,...])
            std::vector<ParticleType> types(n_types);
            std::vector<ParticleID> ids(n_ids);

            for (size_t i = 0; i < n_types; ++i) types[i] = i;
            for (size_t i = 0; i < n_ids; ++i) ids[i] = i;

            // first we gather all interactions (type interactions and id interactions)
            std::vector<InteractionVariant> all_interactions;
            all_interactions.reserve(type_interactions.size() + id_interactions.size());

            for (const auto & interaction : type_interactions) all_interactions.push_back(interaction);
            for (const auto & interaction : id_interactions) all_interactions.push_back(interaction);

            // and we also create a corresponding properties vector for every interaction
            std::vector<InteractionDescriptor> all_interaction_props;
            all_interaction_props.reserve(type_interactions.size() + id_interactions.size());
            for (const auto & interaction : all_interactions) all_interaction_props.push_back(get_properties(interaction));

            // loop through all possible type pairs and register them in the properties of their selected interaction
            for (ParticleType i = 0; i < static_cast<ParticleType>(n_types); i++) {
                for (ParticleType j = 0; j < static_cast<ParticleType>(n_types); j++) {
                    all_interaction_props[type_index(i, j)].used_by_types.emplace_back(i, j);
                }
            }

            // loop through all possible (relevant) id pairs and register them in the properties of their selected interaction
            for (ParticleID i = 0; i < static_cast<ParticleID>(n_ids); i++) {
                for (ParticleID j = i+1; j < static_cast<ParticleID>(n_ids); j++) {
                    all_interaction_props[type_interactions.size() + id_index(i, j)].used_by_ids.emplace_back(i, j);
                }
            }

            // now we merge the properties of all identical interactions
            // first we create a vector of unique interactions and track which entry maps to an entry in unique_interactions
            std::vector<size_t> remapping(all_interactions.size());
            std::vector<InteractionVariant> unique_interactions;
            std::vector<InteractionDescriptor> unique_props;

            for (size_t i = 0; i < all_interactions.size(); i++) {
                const auto & current_interaction = all_interactions[i];
                auto& current_prop  = all_interaction_props[i];

                // check if current_interaction is already contained in unique interactions
                bool found = false;
                size_t found_idx = 0;

                for (size_t j = 0; j < unique_interactions.size(); ++j) {
                    if (is_equal(current_interaction, unique_interactions[j])) {
                        found = true;
                        found_idx = j;
                        break;
                    }
                }

                if (found) {
                    // current interaction is a duplicate of unique_interactions[found_idx] -> merge
                    auto & props = unique_props[found_idx];

                    props.used_by_types.insert(props.used_by_types.end(), current_prop.used_by_types.begin(), current_prop.used_by_types.end());
                    props.used_by_ids.insert(props.used_by_ids.end(), current_prop.used_by_ids.begin(), current_prop.used_by_ids.end());

                    remapping[i] = found_idx;
                } else {
                    // current interaction is not in unique_interactions -> create new entry
                    const size_t new_idx = unique_interactions.size();
                    unique_interactions.push_back(current_interaction);
                    unique_props.push_back(std::move(all_interaction_props[i]));

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

            return InteractionMap {
                .types = types,
                .ids = ids,
                .interactions = unique_props,
                .type_interaction_matrix = type_interaction_matrix,
                .id_interaction_matrix = id_interaction_matrix,
            };
        }


        template<typename Func>
        void dispatch(const ParticleType t1, const ParticleType t2, Func && func) const {
            const auto & variant = get_type_interaction(t1, t2);
            std::visit([&]<IsInteraction I>(const I & interaction) -> void {
                if constexpr (!std::same_as<I, InteractionSentinel> && !std::same_as<I, NoInteraction>) {
                    func(interaction);
                }
            }, variant);
        }

        template<typename Func>
        void dispatch_id(const ParticleID id1, const ParticleID id2, Func && func) const {
            const auto & variant = get_id_interaction(id1, id2);
            std::visit([&]<IsInteraction I>(const I & interaction) -> void {
                if constexpr (!std::same_as<I, InteractionSentinel> && !std::same_as<I, NoInteraction>) {
                    func(interaction);
                }
            }, variant);
        }


        [[nodiscard]] bool has_id_interaction(const ParticleID a, const ParticleID b) const noexcept{
            return a < n_ids && b < n_ids;
        }

        InteractionVariant & get_type_interaction(const ParticleType a, const ParticleType b) noexcept {
            return type_interactions[type_index(a, b)];
        }

        InteractionVariant & get_id_interaction(const ParticleID a, const ParticleID b) noexcept {
            return id_interactions[id_index(a, b)];
        }

        const InteractionVariant& get_type_interaction(const ParticleType a, const ParticleType b) const noexcept {
            return type_interactions[type_index(a,b)];
        }

        const InteractionVariant& get_id_interaction(const ParticleID a, const ParticleID b) const noexcept {
            return id_interactions[id_index(a,b)];
        }


    private:
        std::vector<InteractionVariant> type_interactions; // Interactions between different particle types (e.g. type A <-> type B)
        std::vector<InteractionVariant> id_interactions; // Interactions between specific particle instances (by ID e.g. id1 <-> id2)
        size_t n_types{};
        size_t n_ids{};

        double max_cutoff = 0;

        [[nodiscard]] size_t type_index(const ParticleType a, const ParticleType b) const noexcept{
            return n_types * a + b;
        }

        [[nodiscard]] size_t id_index(const ParticleID a, const ParticleID b) const noexcept{
            return n_ids * a + b;
        }



        void build_type_interactions(std::vector<Type_Interaction>& type_infos, const TypeMap & type_map)
        {
            // collect unique particle types to define types map size (implementation types are dense [0, N-1])
            std::unordered_set<ParticleType> particle_types;
            for (auto& x : type_infos) {
                particle_types.insert(type_map.at(x.type1));
                particle_types.insert(type_map.at(x.type2));
            }

            n_types = particle_types.size();
            type_interactions.resize(n_types * n_types);

            // insert type interactions into map & apply user mappings
            for (auto& x : type_infos) {
                const auto a = type_map.at(x.type1);
                const auto b = type_map.at(x.type2);
                type_interactions[type_index(a, b)] = x.interaction;
                type_interactions[type_index(b, a)] = x.interaction;
            }

            //  mix missing type pairs from diagonals
            for (size_t a = 0; a < n_types; ++a) {
                for (size_t b = 0; b < n_types; ++b) {
                    auto& interaction = get_type_interaction(a, b);
                    if (a == b || !std::holds_alternative<InteractionSentinel>(interaction)) continue;

                    auto& interaction_a = get_type_interaction(a, a);
                    auto& interaction_b = get_type_interaction(b, b);

                    auto mixed_interaction = std::visit([]<typename I1, typename I2>(I1 const& a, I2 const& b) -> InteractionVariant {
                            if constexpr (std::same_as<I1, I2>)
                                return a.mix_interactions(b);
                            else
                                throw std::invalid_argument("Cannot mix different interaction types");
                        },
                        interaction_a, interaction_b);

                    type_interactions[type_index(a, b)] = mixed_interaction;
                    type_interactions[type_index(b, a)] = mixed_interaction;
                }
            }
        }


        void build_id_interactions(std::vector<Id_Interaction>& id_infos, const IdMap & id_map)
        {
            // collect particle ids to define ids map size (implementation ids are dense [0, M-1])
            std::unordered_set<ParticleID> ids;
            for (auto& x : id_infos) {
                ids.insert(id_map.at(x.id1));
                ids.insert(id_map.at(x.id2));
            }

            n_ids = ids.size();
            id_interactions.resize(n_ids * n_ids);

            // insert id interactions into map & apply usr mappings
            for (auto& x : id_infos) {
                const auto a = id_map.at(x.id1);
                const auto b = id_map.at(x.id2);
                id_interactions[id_index(a, b)] = x.interaction;
                id_interactions[id_index(b, a)] = x.interaction;
            }

            // Fill undefined id interactions with no interactions
            for (size_t a = 0; a < n_ids; a++) {
                for (size_t b = 0; b < n_ids; b++) {
                    auto & v = id_interactions[id_index(a, b)];
                    if (a != b && std::holds_alternative<InteractionSentinel>(v)) {
                        v = NoInteraction();
                    }
                }
            }
        }

        void validate_interaction_table() const {
            #ifndef NDEBUG
            for (size_t i = 0; i < n_types; ++i)
                for (size_t j = 0; j < n_types; ++j)
                    APRIL_ASSERT(!std::holds_alternative<InteractionSentinel>(type_interactions[type_index(i, j)]),
                              "inter_type_interactions should not contain InteractionSentinel");

            for (size_t i = 0; i < n_ids; ++i)
                for (size_t j = 0; j < n_ids; ++j) {
                    auto& v = id_interactions[id_index(i, j)];
                    if (i == j)
                        APRIL_ASSERT(std::holds_alternative<InteractionSentinel>(v),
                                  "intra_particle_interactions should contain InteractionSentinel for identical ids");
                    else
                        APRIL_ASSERT(!std::holds_alternative<InteractionSentinel>(v),
                                  "intra_particle_interactions should not contain InteractionSentinel for differing ids");
            }
            #endif
        }

    };
} // namespace april::interaction::internal


















