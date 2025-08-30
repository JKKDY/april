#pragma once
#include <cassert>
#include <functional>
#include <utility>
#include <unordered_set>

#include "april/env/particle.h"
#include "april/env/force.h"
#include "april/common.h"
#include "april/env/environment.h"


namespace april::env::impl {

    // internal placeholder only
    struct NullForce {
        double cutoff_radius = -1.0;
        vec3 operator()(const Particle&, const Particle&, const vec3&) const noexcept {
            assert(false && "NullForce should never be executed");
             return {};
        }
        [[nodiscard]] NullForce mix(NullForce const&) const { return {}; }
    };


    template<class Env> class InteractionManager;

    template<IsForce... Fs>
    class InteractionManager<Environment<ForcePack<Fs...>>> {

    public:
        using force_variant_t = std::variant<NullForce, Fs..., NoForce>;
        using force_variant_info_t = std::variant<Fs...>;
        // using info_t = InteractionInfo<force_variant_t>;

		InteractionManager() = default;
        void build(std::vector<InteractionInfo<force_variant_info_t>> & interaction_infos,
            const std::unordered_map<env::ParticleType, impl::ParticleType> & usr_types_to_impl_types,
            const std::unordered_map<env::ParticleID, impl::ParticleID> & usr_ids_to_impl_ids
        );

        [[nodiscard]] vec3 evaluate(const Particle& p1, const Particle& p2) const;
		[[nodiscard]] vec3 evaluate(const Particle& p1, const Particle& p2, const vec3& r) const;

        [[nodiscard]] double get_max_cutoff() const;

    private:
        std::vector<force_variant_t> inter_type_forces; // Forces between different particle types (e.g. type A <-> type B)
        std::vector<force_variant_t> intra_particle_forces; // Forces between specific particle instances (by ID e.g. id1 <-> id2)

        [[nodiscard]] size_t type_index(const size_t a, const size_t b) const noexcept{
            return n_types * a + b;
        }

        [[nodiscard]] size_t id_index(const size_t a, const size_t b) const noexcept{
            return n_ids * a + b;
        }

        force_variant_t & get_type_force(const size_t a, const size_t b) noexcept {
            return inter_type_forces[type_index(a, b)];
        }

        force_variant_t & get_id_force(const size_t a, const size_t b) noexcept {
            return intra_particle_forces[id_index(a, b)];
        }

        const force_variant_t& get_type_force(size_t a, size_t b) const noexcept {
            return inter_type_forces[type_index(a,b)];
        }
        const force_variant_t& get_id_force(size_t a, size_t b) const noexcept {
            return intra_particle_forces[id_index(a,b)];
        }

        size_t n_types{};
        size_t n_ids{};

        double max_cutoff = 0;
    };


    template <IsForce ... Fs>
    void InteractionManager<Environment<ForcePack<Fs...>>>::build(
        std::vector<InteractionInfo<force_variant_info_t>>& interaction_infos,
        const std::unordered_map<env::ParticleType, impl::ParticleType>& usr_types_to_impl_types,
        const std::unordered_map<env::ParticleID, impl::ParticleID>& usr_ids_to_impl_ids) {

        // partition type vs id
        const auto it = std::partition(interaction_infos.begin(), interaction_infos.end(),
           [](const auto& info) {return info.pair_contains_types;});

        // contains all interaction infos for particle types
        std::vector<InteractionInfo<force_variant_info_t>> type_infos{
            std::make_move_iterator(interaction_infos.begin()),
            std::make_move_iterator(it)
        };

        // contains all interaction infos for particle (id) pairs
        std::vector<InteractionInfo<force_variant_info_t>> id_infos{
            std::make_move_iterator(it),
            std::make_move_iterator(interaction_infos.end())
        };


        // helper lambda (local) to widen info-variant -> manager-variant
        auto widen = [](auto&& small) -> force_variant_t {
            using T = std::decay_t<decltype(small)>;
            return force_variant_t{std::forward<T>(small)};
        };


        // collect unique particle types to define types map size (implementation types are dense [0, N-1])
        std::unordered_set<ParticleType> particle_types_set;
        for (auto & x : type_infos) {
            particle_types_set.insert(usr_types_to_impl_types.at(x.key_pair.first));
            particle_types_set.insert(usr_types_to_impl_types.at(x.key_pair.second));
        }
        this->n_types = particle_types_set.size();
        inter_type_forces = std::vector<force_variant_t>(n_types*n_types);

        // collect particle ids to define ids map size (implementation ids are dense [0, M-1])
        std::unordered_set<ParticleID> particle_id_set;
        for (auto & x : id_infos) {
            particle_id_set.insert(usr_ids_to_impl_ids.at(x.key_pair.first));
            particle_id_set.insert(usr_ids_to_impl_ids.at(x.key_pair.second));
        }
        this->n_ids = particle_id_set.size();
        intra_particle_forces = std::vector<force_variant_t>(n_ids*n_ids);


        // insert type forces into map & apply usr mappings
        for (auto & x : type_infos) {
            const ParticleType a = usr_types_to_impl_types.at(x.key_pair.first);
            const ParticleType b = usr_types_to_impl_types.at(x.key_pair.second);

            inter_type_forces[type_index(a, b)] = std::visit(widen, x.force);
            inter_type_forces[type_index(b, a)] = std::visit(widen, x.force);
        }

        //  mix missing type pairs from diagonals
        for (size_t a = 0; a < n_types; a++) {
            for (size_t b = 0; b < n_types; b++) {
                auto & force = get_type_force(a, b);
                if (a == b || !std::holds_alternative<NullForce>(force)) continue;

                auto &va = get_type_force(a,a);
                auto &vb = get_type_force(b,b);

                auto f = std::visit([]<typename T0, typename T1>(T0 const& fa, T1 const& fb)->force_variant_t {
                    using A = std::decay_t<T0>;
                    using B = std::decay_t<T1>;
                    if constexpr (std::same_as<A,B>) {
                        return fa.mix(fb); // returns A
                    } else {
                        throw std::invalid_argument("Cannot mix different force types");
                    }
                }, va, vb);

                inter_type_forces[type_index(a, b)] = f;
                inter_type_forces[type_index(b, a)] = f;
            }
        }


        // insert id forces into map & apply usr mappings
        for (auto & x : id_infos) {
            const ParticleType a = usr_ids_to_impl_ids.at(x.key_pair.first);
            const ParticleType b = usr_ids_to_impl_ids.at(x.key_pair.second);

            intra_particle_forces[id_index(a, b)] = std::visit(widen, x.force);
            intra_particle_forces[id_index(b, a)] = std::visit(widen, x.force);
        }

        // Fill undefined id interactions with no forces
        for (size_t a = 0; a < n_ids; a++) {
            for (size_t b = 0; b < n_ids; b++) {
                auto & v = intra_particle_forces[id_index(a, b)];
                if (a != b && std::holds_alternative<NullForce>(v)) {
                    v = NoForce();
                }
            }
        }


        // check if force maps are valid
        for (size_t i = 0; i < n_types; i++) {
            for (size_t j = 0; j < n_types; j++) {
                AP_ASSERT(!std::holds_alternative<NullForce>(inter_type_forces[type_index(i, j)]),
                    "inter_type_forces should not contain NullForce");
            }
        }

        for (size_t i = 0; i < n_ids; i++) {
            for (size_t j = 0; j < n_ids; j++) {
                if (i == j) {
                    AP_ASSERT(std::holds_alternative<NullForce>(intra_particle_forces[id_index(i, j)]),
                        "intra_particle_forces should contain NullForce for p1.id = p2.id");
                } else {
                    AP_ASSERT(!std::holds_alternative<NullForce>(intra_particle_forces[id_index(i, j)]),
                        "intra_particle_forces should not contain NullForce for differing particle ids");
                }
            }
        }


        // get the max cutoff distance
        auto cutoff_of = [](auto const& v){
            return std::visit([](auto const& f){ return f.cutoff_radius; }, v);
        };

        max_cutoff = 0.0;
        for (auto const& v : inter_type_forces)
            max_cutoff = std::max(max_cutoff, cutoff_of(v));

        for (auto const& v : intra_particle_forces)
            if (!std::holds_alternative<NullForce>(v))
                max_cutoff = std::max(max_cutoff, cutoff_of(v));
    }

    template <IsForce ... Fs>
    vec3 InteractionManager<Environment<ForcePack<Fs...>>>::evaluate(
        const Particle& p1, const Particle& p2) const {
        return evaluate(p1, p2, p2.position - p1.position); // dist vector points from p1 to p2
    }

    template <IsForce ... Fs>
    vec3 InteractionManager<Environment<ForcePack<Fs...>>>::evaluate(
        const Particle& p1, const Particle& p2, const vec3& r) const {

        auto & tF = get_type_force(p1.type, p2.type);
        vec3 force = std::visit([&](auto const& f){ return f(p1,p2,r); }, tF);

        // check if both particles even have any individual interactions defined for them
        if (p1.id < n_ids && p2.id < n_ids) {
             auto & iF = get_id_force(p1.id, p2.id);
             force += std::visit([&](auto const& f){ return f(p1,p2,r); }, iF);
        }

        return force;
    }

    template <IsForce ... Fs>
    double InteractionManager<Environment<ForcePack<Fs...>>>::get_max_cutoff() const {
        return max_cutoff;
    }


} // namespace april::env::impl



