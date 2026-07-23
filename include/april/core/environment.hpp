/**
* @file environment.hpp
 * @brief Defines APRIL's declarative simulation-configuration interface.
 *
 * Environment collects particles, interactions, boundaries, controllers,
 * fields, particle attributes, and domain settings before runtime resources
 * are created. It acts as a staging area.
 *
 * An Environment does not own a running simulation. Its contents are validated
 * and materialized into a System by build_system().
 */

#pragma once

#include <any>
#include <array>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "april/base/types.hpp"
#include "april/core/domain.hpp"
#include "april/core/internal/environment_traits.hpp"

#include "april/interactions/force.hpp"
#include "april/boundaries/boundary.hpp"
#include "april/controllers/controller.hpp"
#include "april/fields/field.hpp"

#include "april/particle/generators.hpp"
#include "april/particle/properties.hpp"
#include "april/particle/attributes.hpp"


namespace april {

    /**
     * @brief Selects interactions within one particle type.
     *
     * When passed to Environment::add_interaction(), the force is registered for
     * pairs in which both particles have `type`.
     */
    struct to_type {
        ParticleType type;  /// Particle type participating in the self-interaction.
        constexpr explicit to_type(const ParticleType type) noexcept
            : type(type) {}
    };

    /**
     * @brief Selects interactions between two particle types.
     *
     * The ordering of `t1` and `t2` does not define an interaction direction unless
     * the configured force explicitly has directional semantics.
     */
    struct between_types {
        ParticleType t1;
        ParticleType t2;

        constexpr between_types(const ParticleType first, const ParticleType second) noexcept
            : t1(first),
              t2(second)
        {}
    };

    /**
     * @brief Selects an interaction between two persistent particle identifiers.
     *
     * ID-scoped interactions are commonly used for bonds and other topology-based
     * relationships.
     */
    struct between_ids {
        ParticleID id1;
        ParticleID id2;

        constexpr between_ids(const ParticleID first, const ParticleID second) noexcept
            : id1(first),
              id2(second)
        {}
    };


    /**
      * @brief Collects a declarative APRIL simulation definition.
      *
      * Environment stores initial particles, interaction instances, boundary
      * conditions, controllers, fields, particle attributes, and domain settings.
      * Calling build_system() validates and transforms this definition into the
      * concrete resources owned by a running System.
      *
      * The component packs specify which concrete component types may be added.
      * They do not themselves add component instances to the environment.
      *
      * @tparam FPack Force types accepted by add_interaction().
      * @tparam BPack Boundary types accepted by set_boundary() and set_boundaries().
      * @tparam CPack Controller types accepted by add_controller().
      * @tparam FFPack Field types accepted by add_field().
      * @tparam ParticleAttributes User-defined per-particle attribute definition.
      */
    template<
        interactions::internal::IsForcePack FPack,
        boundary::internal::IsBoundaryPack BPack,
        controller::internal::IsControllerPack CPack,
        field::internal::IsFieldPack FFPack,
        particle::IsParticleAttributes ParticleAttributes>
    class Environment {
    public:
        using traits = core::internal::EnvironmentTraits<
            FPack,
            BPack,
            CPack,
            FFPack,
            ParticleAttributes
        >;

        /**
         * @brief Constructs an empty environment with the specified component packs.
         *
         * The constructor arguments are type markers used to select the environment's
         * supported component and particle-attribute types. Their values are not
         * retained.
         */
        explicit Environment(
            FPack,
            BPack,
            CPack,
            FFPack,
            ParticleAttributes
        ) {}

        /// Constructs an empty environment with no custom component types or attributes.
        Environment()
            : Environment(
                forces<>,
                boundaries<>,
                controllers<>,
                fields<>,
                NoParticleAttributes{}
            )
        {}

        /**
         * @brief Constructs an empty environment and deduces its supported types.
         *
         * @code
         * auto environment = Environment(
         *     forces<LennardJones>,
         *     boundaries<ReflectiveBoundary>
         * );
         * @endcode
         *
         * Constructor arguments must be APRIL component packs or a particle-attribute
         * definition. Missing pack categories are replaced by empty packs.
         *
         * @tparam Args Component-pack and particle-attribute marker types.
         * @param args Type-selection markers. Their runtime values are not retained.
         */
        template<class... Args>
        requires (core::internal::is_any_pack_v<std::remove_cvref_t<Args>> && ...) &&
            (!std::same_as<std::remove_cvref_t<Args>, Environment> && ...)
        explicit Environment(Args&&...)
            : Environment(
                core::internal::get_pack_t<interactions::internal::ForcePack, Args...>{},
                core::internal::get_pack_t<boundary::internal::BoundaryPack, Args...>{},
                core::internal::get_pack_t<controller::internal::ControllerPack,Args...>{},
                core::internal::get_pack_t<field::internal::FieldPack, Args...>{},
                core::internal::get_particle_attributes_t<Args...>{}
            ) {}

    private:
        traits::environment_data_t data;
        friend auto core::internal::get_env_data<FPack, BPack, CPack, FFPack> (const Environment& env);

    public:

        //--------------
        // ADD PARTICLES
        //--------------
        /**
         * @brief Adds one particle to the environment.
         *
         * A user-specified identifier must be unique among all particles already
         * present in the environment. Particles without identifiers receive IDs during
         * the system build process.
         *
         * @param particle Particle definition to add.
         *
         * @throws std::invalid_argument If `particle.id` is present and already used.
         */
        void add_particle(const Particle& particle) {
            if (particle.id.has_value() && data.user_particle_ids.contains(particle.id.value())) {
                throw std::invalid_argument("specified id is not unique");
            }
            data.particles.push_back(particle);

            data.user_particle_types.insert(particle.type);
            if (particle.id.has_value()) {
                data.user_particle_ids.insert(particle.id.value());
            }
        }

        /**
         * @brief Constructs and adds one active particle.
         *
         * @param position Initial particle position.
         * @param velocity Initial particle velocity.
         * @param mass Initial particle mass.
         * @param type Particle type identifier.
         * @param id Optional persistent particle identifier.
         * @param user_data Optional type-erased user data consumed while building the
         * configured particle attributes.
         *
         * @throws std::invalid_argument If `id` is present and already used.
         */
        void add_particle(
            const vec3& position,
            const vec3& velocity,
            const double mass,
            const ParticleType type = 0,
            const std::optional<ParticleID> id = {},
            std::any user_data = {}
        ) {
            Particle p;
            p.id = id;
            p.type = type;
            p.position = position;
            p.velocity = velocity;
            p.mass = mass;
            p.state = ParticleState::ALIVE;
            p.user_data = std::move(user_data);
            add_particle(p);
        }

        /**
         * @brief Adds a collection of particle definitions.
         *
         * Particles are validated and inserted in input order.
         *
         * @param particles Particle definitions to add.
         *
         * @throws std::invalid_argument If any explicitly specified identifier is
         * duplicated.
         *
         * @note If insertion fails partway through, particles inserted before the
         * failing particle remain in the environment.
         */
        void add_particles(const std::vector<Particle>& particles) {
            this->data.particles.reserve(this->data.particles.size() + particles.size());
            for (const auto& particle : particles)
                add_particle(particle);
        }

        /**
         * @brief Generates and adds particles from a particle generator.
         *
         * @tparam G Particle-generator type.
         * @param generator Generator whose output is added to the environment.
         *
         * @throws std::invalid_argument If generated particles contain duplicate
         * explicitly specified identifiers.
         */
        template<IsParticleGenerator G>
        void add_particles(const G& generator) {
            add_particles(generator.to_particles());
        }


        //-----------------
        // ADD INTERACTIONS
        //-----------------
        /**
         * @brief Registers a force for pairs of one particle type.
         *
         * @tparam F Force type declared in this environment's force pack.
         * @param force Force instance to store.
         * @param scope Particle type to which the self-interaction applies.
         */
        template<interactions::IsForce F>
        requires traits::template is_valid_force_v<F>
        void add_interaction(F && force, to_type scope) {
            data.type_interactions.emplace_back(scope.type, scope.type, typename traits::force_variant_t{std::move(force)});
        }

        /**
         * @brief Registers a force between two particle types.
         *
         * @tparam F Force type declared in this environment's force pack.
         * @param force Force instance to store.
         * @param scope Pair of particle types to which the force applies.
         */
        template<interactions::IsForce F> requires traits::template is_valid_force_v<F>
        void add_interaction(F force, between_types scope) {
            data.type_interactions.emplace_back(scope.t1, scope.t2, typename traits::force_variant_t{std::move(force)});
        }

        /**
         * @brief Registers a force between two specific particles.
         *
         * ID-scoped interactions are stored as topology interactions and remain tied to
         * persistent identifiers rather than physical particle indices.
         *
         * @tparam F Force type declared in this environment's force pack.
         * @param force Force instance to store.
         * @param scope Persistent particle identifiers defining the interaction.
         */
        template<interactions::IsForce F> requires traits::template is_valid_force_v<F>
        void add_interaction(F force, between_ids scope) {
            data.id_interactions.emplace_back(scope.id1, scope.id2, typename traits::force_variant_t{std::move(force)});
        }


        //---------------
        // ADD BOUNDARIES
        //---------------
        /**
         * @brief Sets the boundary condition for one domain face.
         *
         * Replaces any boundary condition currently assigned to the selected face.
         *
         * @tparam B Boundary type declared in this environment's boundary pack.
         * @param boundary Boundary instance to store.
         * @param face Domain face to configure.
         */
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundary(B boundary, const DomainFace face) {
            data.boundaries[boundary::face_to_int(face)].template emplace<B>(std::move(boundary));
        }

        /**
         * @brief Sets the same boundary condition on multiple domain faces.
         *
         * A copy of `boundary` is stored for each selected face. Any existing boundary
         * on those faces is replaced.
         *
         * @tparam B Copy-constructible boundary type declared in this environment's
         * boundary pack.
         * @param boundary Boundary instance copied to each face.
         * @param faces Domain faces to configure.
         */
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundaries(B boundary, const std::vector<DomainFace> & faces) {
            for (const DomainFace face : faces) {
                data.boundaries[boundary::face_to_int(face)].template emplace<B>(boundary);
            }
        }

        /**
         * @brief Sets one boundary condition for each domain face.
         *
         * Array order is `-X`, `+X`, `-Y`, `+Y`, `-Z`, `+Z`.
         *
         * @tparam B Boundary type declared in this environment's boundary pack.
         * @param boundaries Boundary instances indexed in domain-face order.
         */
        template<boundary::IsBoundary B> requires traits::template is_valid_boundary_v<B>
        void set_boundaries(const std::array<B, 6> & boundaries) {
            for (const DomainFace face : all_faces) {
                data.boundaries[boundary::face_to_int(face)].template emplace<B>(boundaries[boundary::face_to_int(face)]);
            }
        }


        //----------------
        // ADD CONTROLLERS
        //----------------
        /**
         * @brief Adds a controller instance.
         *
         * @tparam C Controller type declared in this environment's controller pack.
         * @param controller Controller instance to store.
         */
        template<controller::IsController C>
        requires traits::template is_valid_controller_v<std::remove_cvref_t<C>>
        void add_controller(C && controller) {
            data.controllers.add(std::forward<C>(controller));
        }

        /**
         * @brief Adds multiple controller instances.
         *
         * @tparam C Controller types declared in this environment's controller pack.
         * @param controllers Controller instances to store, in argument order.
         */
        template<controller::IsController... C>
        requires (traits::template is_valid_controller_v<std::remove_cvref_t<C>> && ...)
        void add_controllers(C&&... controllers) {
            (data.controllers.add(std::forward<C>(controllers)), ...);
        }



        //-----------
        // ADD FIELDS
        //-----------
        /**
         * @brief Adds a global field instance.
         *
         * @tparam F Field type declared in this environment's field pack.
         * @param field Field instance to store.
         */
        template<field::IsField F>
        requires traits::template is_valid_field_v<std::remove_cvref_t<F>>
        void add_field(F && field) {
            data.fields.add(std::forward<F>(field));
        }

        /**
         * @brief Adds multiple field instances.
         *
         * @tparam F Field types declared in this environment's field pack.
         * @param fields Field instances to store, in argument order.
         */
        template<field::IsField... F>
        requires (traits::template is_valid_field_v<std::remove_cvref_t<F>> && ...)
        void add_fields(F&&... fields) {
            (data.fields.add(std::forward<F>(fields)), ...);
        }


        //-----------
        // SET DOMAIN
        //-----------
        /**
         * @brief Sets the minimum corner of the simulation domain.
         *
         * @param origin Coordinates of the domain's minimum corner.
         */
        void set_origin(const vec3d& origin) {
            this->data.domain.origin = origin;
        }
        /// @brief Set the coordinate of the minimum corner from scalars.
        void set_origin(const double x, const double y, const double z) {
            set_origin({x,y,z});
        }

        /**
         * @brief Sets the dimensions of the simulation domain.
         *
         * @param extent Domain width along each coordinate axis.
         *
         * @pre Every component of `extent` must be non-negative.
         */
        void set_extent(const vec3d& extent) {
            this->data.domain.extent = extent;
        }
        /// @brief Set the dimensions of the simulation box from scalars.
        void set_extent(const double x, const double y, const double z) {
            set_extent({x,y,z});
        }

        /**
         * @brief Replaces the complete simulation-domain definition.
         *
         * @param domain Domain geometry to store.
         */
        void set_domain(const Domain& domain) {
            data.domain = domain;
        }

        /**
         * @brief Sets absolute padding used when deriving a domain from the particles.
         *
         * During build, the particle bounding box is expanded by the specified amount
         * on each side of each axis.
         *
         * @param padding Absolute per-axis padding.
         */
        void set_domain_padding(const vec3d& padding) {
            data.margin_abs = padding;
        }
        /// Sets the same absolute domain padding on every axis.
        void set_domain_padding(const double padding) {
            data.margin_abs = vec3d{padding};
        }

        /**
         * @brief Sets relative padding used when deriving the domain from particles.
         *
         * Each component is multiplied by the corresponding extent of the particle
         * bounding box and added to both sides of that axis.
         *
         * For example, a factor of `0.5` adds half of the particle extent to each side,
         * producing a final domain extent equal to twice the particle extent.
         *
         * When both absolute and relative padding are configured, the larger resulting
         * padding is selected independently for each axis.
         *
         * @param factor Relative per-axis padding factor.
         */
        void set_domain_padding_factor(const vec3d& factor) {
            data.margin_fac = factor;
        }
        /// Sets the same relative domain-padding factor on every axis.
        void set_domain_padding_factor(const double factor) {
            set_domain_padding_factor(vec3d{factor});
        }


        // ------------------
        // CHAINING INTERFACE
        // ------------------
        //
        // The with_* functions are value-category-preserving wrappers around the
        // corresponding add_*, set_*, and domain-configuration functions.

        /**
         * @name Chaining interface
         *
         * Value-category-preserving alternatives to the corresponding mutating
         * operations. Each function applies the requested change and returns the same
         * Environment object, allowing configuration calls to be chained.
         *
         * @warning A reference returned from chaining on a temporary must not be stored
         * beyond the full expression.
         *
         * @{
         */
        auto&& with_particle(this auto&& self, const Particle& p) {
            self.add_particle(p);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_particle(
            this auto&& self,
            const vec3& position,
            const vec3& velocity,
            double mass,
            ParticleType type = 0,
            std::optional<ParticleID> id = {},
            std::any user_data = {}
        ) {
            self.add_particle(position, velocity, mass, type, id, std::move(user_data));
            return std::forward<decltype(self)>(self);
        }

        auto&& with_particles(this auto&& self, const std::vector<Particle>& ps) {
            self.add_particles(ps);
            return std::forward<decltype(self)>(self);
        }

        template<IsParticleGenerator G>
        auto&& with_particles(this auto&& self, const G& particles) {
            self.add_particles(particles);
            return std::forward<decltype(self)>(self);
        }

        template<interactions::IsForce F>
            requires traits::template is_valid_force_v<std::remove_cvref_t<F>>
        auto&& with_interaction(this auto&& self, F&& force, to_type scope) {
            self.add_interaction(std::forward<F>(force), scope);
            return std::forward<decltype(self)>(self);
        }

        template<interactions::IsForce F>
            requires traits::template is_valid_force_v<std::remove_cvref_t<F>>
        auto&& with_interaction(this auto&& self, F&& force, between_types scope) {
            self.add_interaction(std::forward<F>(force), scope);
            return std::forward<decltype(self)>(self);
        }

        template<interactions::IsForce F>
            requires traits::template is_valid_force_v<std::remove_cvref_t<F>>
        auto&& with_interaction(this auto&& self, F&& force, between_ids scope) {
            self.add_interaction(std::forward<F>(force), scope);
            return std::forward<decltype(self)>(self);
        }

        template<boundary::IsBoundary B>
            requires traits::template is_valid_boundary_v<std::remove_cvref_t<B>>
        auto&& with_boundary(this auto&& self, B&& boundary, DomainFace face) {
            self.set_boundary(std::forward<B>(boundary), face);
            return std::forward<decltype(self)>(self);
        }

        template<boundary::IsBoundary B>
            requires traits::template is_valid_boundary_v<std::remove_cvref_t<B>>
        auto&& with_boundaries(this auto&& self, B&& boundary, const std::vector<DomainFace>& faces) {
            self.set_boundaries(std::forward<B>(boundary), faces);
            return std::forward<decltype(self)>(self);
        }

        template<boundary::IsBoundary B>
            requires traits::template is_valid_boundary_v<B>
        auto&& with_boundaries(this auto&& self, const std::array<B, 6>& boundaries) {
            self.set_boundaries(boundaries);
            return std::forward<decltype(self)>(self);
        }

        template<controller::IsController C>
            requires traits::template is_valid_controller_v<std::remove_cvref_t<C>>
        auto&& with_controller(this auto&& self, C&& controller) {
            self.add_controller(std::forward<C>(controller));
            return std::forward<decltype(self)>(self);
        }

        template<controller::IsController... C>
        requires (traits::template is_valid_controller_v<std::remove_cvref_t<C>> && ...)
        auto&& with_controllers(this auto&& self, C&&... controllers) {
            self.add_controllers(std::forward<C>(controllers)...);
            return std::forward<decltype(self)>(self);
        }

        template<field::IsField F>
        requires traits::template is_valid_field_v<std::remove_cvref_t<F>>
        auto&& with_field(this auto&& self, F&& field) {
            self.add_field(std::forward<F>(field));
            return std::forward<decltype(self)>(self);
        }

        template<field::IsField... F>
        requires (traits::template is_valid_field_v<std::remove_cvref_t<F>> && ...)
        auto&& with_fields(this auto&& self, F&&... fields) {
            self.add_fields(std::forward<F>(fields)...);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_origin(this auto&& self, const vec3& o) {
            self.set_origin(o);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_extent(this auto&& self, const vec3& e) {
            self.set_extent(e);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_origin(this auto&& self, double x, double y, double z) {
            self.set_origin(x, y, z);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_extent(this auto&& self, double x, double y, double z) {
            self.set_extent(x, y, z);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_domain(this auto&& self, const Domain& domain) {
            self.set_domain(domain);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_domain_padding(this auto&& self, double margin) {
            self.set_domain_padding(margin);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_domain_padding(this auto&& self, const vec3& margin) {
            self.set_domain_padding(margin);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_domain_padding_factor(this auto&& self, double factor) {
            self.set_domain_padding_factor(factor);
            return std::forward<decltype(self)>(self);
        }

        auto&& with_domain_padding_factor(this auto&& self, const vec3& factor) {
            self.set_domain_padding_factor(factor);
            return std::forward<decltype(self)>(self);
        }

        /** @} */
    };


    /**
     * @brief CTAD Guide to deduce template parameters from constructor arguments.
     * Allows: auto env = Environment(forces<LJ>, boundaries<Reflective>);
     */
    template<class... Args>
    Environment(Args...)
        -> Environment<
            core::internal::get_pack_t<interactions::internal::ForcePack, Args...>,
            core::internal::get_pack_t<boundary::internal::BoundaryPack, Args...>,
            core::internal::get_pack_t<controller::internal::ControllerPack,Args...>,
            core::internal::get_pack_t<field::internal::FieldPack, Args...>,
            core::internal::get_particle_attributes_t<Args...>
        >;


    namespace core {
        namespace internal {
            template<typename T>
            inline constexpr bool is_environment_v = false;


            // Recognize concrete Environment specializations.
            template<
                interactions::internal::IsForcePack FPack,
                boundary::internal::IsBoundaryPack BPack,
                controller::internal::IsControllerPack CPack,
                field::internal::IsFieldPack FFPack,
                particle::IsParticleAttributes ParticleAttributes>
            inline constexpr bool is_environment_v<
                Environment<
                    FPack,
                    BPack,
                    CPack,
                    FFPack,
                    ParticleAttributes
                >
            > = true;
        }

        /**
         * @brief Identifies concrete APRIL Environment types.
         *
         * @tparam T Type to test.
         */
        template<typename T>
        concept IsEnvironment = internal::is_environment_v<std::remove_cvref_t<T>>;
    }
}



