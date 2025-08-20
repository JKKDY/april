#include "april/env/environment.h"
#include <algorithm>
#include <stdexcept>
#include <ranges>
#include <unordered_set>

#include "april/april.h"


namespace april::env {
    void Environment::add_particle(const vec3& position, const vec3& velocity, const double mass, const ParticleType type, const ParticleID id) {
        add_particle(Particle{
            .id = id,
            .type = type,
            .position =  position,
            .velocity = velocity,
            .mass = mass,
            .state = ParticleState::ALIVE
        });
    }

    void Environment::add_particle(const Particle & particle) {
        if (particle.id != PARTICLE_ID_DONT_CARE && data.usr_particle_ids.contains(particle.id)) {
            throw std::invalid_argument("specified id is not unique");
        }

        data.particles.push_back(particle);

        if (particle.id != PARTICLE_ID_DONT_CARE) {
            data.usr_particle_ids.insert(particle.id);
        }

        data.usr_particle_types.insert(particle.type);
    }

    void Environment::add_particles(const std::vector<Particle> & particles) {
        this->data.particles.reserve(this->data.particles.size() + particles.size());

        for (auto & p : particles) {
            add_particle(p);
        }
    }

     std::vector<ParticleID> Environment::add_particle_cuboid(const ParticleCuboid& cuboid) {
        const uint32_t particle_count = cuboid.particle_count[0] * cuboid.particle_count[1] * cuboid.particle_count[2];
        const double width = cuboid.distance;

        std::vector<ParticleID> ids;
        ids.reserve(particle_count);

        const auto it = std::ranges::max_element(
            data.particles,
            {},               // default `<` comparator
            &Particle::id       // project each Particle to its `id`
        );
        int id = ( it == data.particles.end() ? 0 : it->id+1 );

        data.particles.reserve(data.particles.size() + particle_count);

        for (unsigned int x = 0; x < cuboid.particle_count[0]; ++x) {
            for (unsigned int y = 0; y < cuboid.particle_count[1]; ++y) {
                for (unsigned int z = 0; z < cuboid.particle_count[2]; ++z) {

                    ids.push_back(id);

                    Particle p = {
                        .id = id++,
                        .type = cuboid.type,
                        .position = cuboid.origin + vec3(x * width, y * width, z * width),
                        .velocity = cuboid.mean_velocity,
                        .mass = cuboid.mass,
                        .state = cuboid.state,
                    };
                    p.velocity += cuboid.thermal_velocity(p);

                    add_particle(p);
                }
            }
        }

        return ids;
    }


    std::vector<ParticleID> Environment::add_particle_sphere(const ParticleSphere& sphere) {
        const double width = sphere.distance;
        const vec3 & r = sphere.radius;

        std::vector<ParticleID> ids;

        // get the maximum current id
        const auto it = std::ranges::max_element(
            data.particles,
            {},                 // default `<` comparator
            &Particle::id       // project each Particle to its `id`
        );
        int id = (it == data.particles.end() ? 0 : it->id+1);

        for (int x = -static_cast<int>(sphere.radius.x/width); x < sphere.radius.x; ++x) {
            for (int y = -static_cast<int>(sphere.radius.y/width); y < sphere.radius.y; ++y) {
                for (int z = -static_cast<int>(sphere.radius.z/width); z < sphere.radius.z; ++z) {

                    vec3 pos = {x * width, y * width, z * width};
                    vec3 pos_sq = pos.mul(pos_sq);

                    // if not in ellipsoid skip
                    if ( pos_sq.x / r.x + pos_sq.y / r.y + pos_sq.z / r.z > 1) continue;

                    ids.push_back(id);
                    Particle p =  {
                        .id = id++,
                        .type = sphere.type,
                        .position = sphere.center + pos,
                        .velocity = sphere.mean_velocity,
                        .mass = sphere.mass,
                        .state = sphere.state,
                    };
                    p.velocity += sphere.thermal_velocity(p);

                    add_particle(p);
                }
            }
        }
        return ids;
    }

    void Environment::set_origin(const vec3 & origin) {
        this->data.domain.origin = origin;
    }
    void Environment::set_extent(const vec3 & extent) {
        if (extent.x < 0 || extent.y < 0 || extent.z < 0) {
            throw std::invalid_argument("Extent cannot be negative. Got: " + extent.to_string());
        }
        this->data.domain.extent = extent;
    }

    void Environment::add_force_field() {
        throw std::logic_error("Not implemented yet");
    }
    void Environment::add_barrier() {
        throw std::logic_error("Not implemented yet");
    }
    void Environment::set_boundary_conditions() {
        throw std::logic_error("Not implemented yet");
    }


    namespace impl {
        EnvironmentData get_env_data(const Environment& env) {
            return env.data;
        }
    }
}



namespace april::env {
    void Environment::update_forces() const {
        container->calculate_forces();
    }

    impl::ParticleIterator Environment::particles(const ParticleState state) {
         return {particle_storage, state};
    }
} // namespace april::env

