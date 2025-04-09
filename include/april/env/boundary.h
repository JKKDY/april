#pragma once

#include <functional>

#include "Common.h"


namespace april::env {
	struct Particle;
	struct GridCell;

    struct BoundaryNormal {
        static const int3 LEFT;    ///< Normal vector of the left face.
        static const int3 RIGHT;   ///< Normal vector of the right face.
        static const int3 TOP;     ///< Normal vector of the top face.
        static const int3 BOTTOM;  ///< Normal vector of the bottom face.
        static const int3 FRONT;   ///< Normal vector of the front face.
        static const int3 BACK;    ///< Normal vector of the back face.
    };

    // Enumeration of possible boundary conditions.
    enum BoundaryRule {
        OUTFLOW,             ///< Particles that cross the boundary are being removed.
        PERIODIC,            ///< Particles exiting one side reenter from the opposite side.
        REPULSIVE_FORCE,     ///< Particles near the boundary are subject to a repulsive force.
        VELOCITY_REFLECTION, ///< When a particle crosses the boundary, its velocity is reflected about the boundary.
    };

    using BoundaryForce = std::function<double(double)>;


    class  Boundary {
    public:
        enum Extent {
            WIDTH,   ///< Width of the simulation space.
            HEIGHT,  ///< Height of the simulation space.
            DEPTH    ///< Depth of the simulation space.
        };

        enum Face {
            LEFT,
            RIGHT,
            TOP,
            BOTTOM,
            FRONT,
            BACK
        };

        static Face normal_to_face(const int3& normal);

      
        Boundary();
        void set_boundary_rule(BoundaryRule rule);
        void set_boundary_rule(BoundaryRule rule, const int3& face_normal);


        void set_boundary_force(const BoundaryForce& force);

        void apply_boundary(Particle& particle, const GridCell& current_cell, const GridCell& previous_cell) const;

        [[nodiscard]] const std::array<BoundaryRule, 6>& boundary_rules() const;


      
        static BoundaryForce LennardJonesForce(double epsilon, double sigma);
        static BoundaryForce InverseDistanceForce(double cutoff, double pre_factor, int exponent = 2);

        // functions for querying whether the force has been set if required.
        [[nodiscard]] bool requires_force_function() const;
        [[nodiscard]] bool has_force_function() const;

        vec3 extent{ MAX_EXTENT, MAX_EXTENT, MAX_EXTENT }; ///< Dimensions of the boundary [width, height, depth].
        vec3 origin{ 0, 0, 0 }; ///< origin of the boundary.
    private:
        void apply_rule(const int3& normal, Particle& particle, const GridCell& current_cell) const;
        void outflow_rule(Particle& particle, const GridCell& previous_cell) const;
        void periodic_rule(Particle& particle, const int3& normal, const GridCell& current_cell) const;
        void repulsive_force_rule(Particle& particle, const int3& normal, const GridCell& cell) const;
        void velocity_reflection_rule(Particle& particle, const int3& normal, const GridCell& cell) const;

        std::array<BoundaryRule, 6> rules{}; ///< Rules applied to each face of the boundary [left, right, top, bottom, front, back]
        BoundaryForce force;  ///< The applied force at the boundary.
    };

}