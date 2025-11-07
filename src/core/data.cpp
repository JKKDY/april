#include "april/env/data.h"
#include "april/env/environment.h"

#include <algorithm>
#include <ranges>

namespace april::env::internal {

	void add_particle_impl(EnvironmentCommonData& data, const Particle& particle) {
		if (particle.id.has_value() && data.user_particle_ids.contains(particle.id.value())) {
			throw std::invalid_argument("specified id is not unique");
		}

		data.particles.push_back(particle);

		if (particle.id.has_value()) {
			data.user_particle_ids.insert(particle.id.value());
		}

		data.user_particle_types.insert(particle.type);
	}


	std::vector<ParticleID> add_cuboid_particles_impl(EnvironmentCommonData& data, const ParticleCuboid& cuboid) {
		if (cuboid.distance == 0) {
			throw std::logic_error("Cuboid inter-particle distance is set to 0!");
		}
		const uint32_t particle_count = cuboid.particle_count[0] * cuboid.particle_count[1] * cuboid.particle_count[2];
		const double width = cuboid.distance;

		std::vector<ParticleID> ids;
		ids.reserve(particle_count);

		auto valid = data.particles | std::views::filter([](auto& p){ return p.id.has_value(); });
		const auto it = std::ranges::max_element (
			valid,
			{},               // default `<` comparator
			&Particle::id       // project each Particle to its `id`
		);
		ParticleID id = ( it == valid.end() ? 0 : it->id.value()+1 );

		data.particles.reserve(data.particles.size() + particle_count);

		for (unsigned int x = 0; x < cuboid.particle_count[0]; ++x) {
			for (unsigned int y = 0; y < cuboid.particle_count[1]; ++y) {
				for (unsigned int z = 0; z < cuboid.particle_count[2]; ++z) {

					ids.push_back(id);

					Particle p;

					p.id = id++;
					p.type		= cuboid.type_idx;
					p.position	= cuboid.origin + vec3(x * width, y * width, z * width);
					p.velocity	= cuboid.mean_velocity;
					p.mass		= cuboid.particle_mass;
					p.state		= cuboid.particle_state;
					p.user_data	= cuboid.user_data;
					p.velocity += cuboid.thermal_velocity(p.position);

					add_particle_impl(data, p);
				}
			}
		}
		return ids;
	}


	std::vector<ParticleID> add_sphere_particles_impl(EnvironmentCommonData& data, const ParticleSphere& sphere) {
		if (sphere.distance == 0) {
            throw std::logic_error("Sphere inter-particle distance is set to 0!");
        }
        std::vector<ParticleID> ids;

        // get the maximum current id
		auto valid = data.particles | std::views::filter([](auto& p){ return p.id.has_value(); });
		const auto it = std::ranges::max_element (
			valid,
			{},               // default `<` comparator
			&Particle::id       // project each Particle to its `id`
		);
		ParticleID id = ( it == valid.end() ? 0 : it->id.value()+1 );

        const double width = sphere.distance;
        const vec3 radii = {
            std::max(sphere.radii.x, width),
            std::max(sphere.radii.y, width),
            std::max(sphere.radii.z, width)
        };

        for (int x = -static_cast<int>(radii.x/width); x < static_cast<int>(radii.x/width); ++x) {
            for (int y = -static_cast<int>(radii.y/width); y < static_cast<int>(radii.y/width); ++y) {
                for (int z = -static_cast<int>(radii.z/width); z < static_cast<int>(radii.z/width); ++z) {

                    const vec3 pos = {x * width, y * width, z * width};
                    const vec3 pos_sq = pos * pos;

                    // if not in ellipsoid skip
                    if (pos_sq.x/(radii.x*radii.x) +
                        pos_sq.y/(radii.y*radii.y) +
                        pos_sq.z/(radii.z*radii.z) >= 1) continue;

                    ids.push_back(id);

                	Particle p;
                    p.id = id++;
                    p.type = sphere.type_idx;
                    p.position = sphere.center + pos;
                    p.velocity = sphere.mean_velocity;
                    p.mass = sphere.particle_mass;
                    p.state = sphere.particle_state;
                	p.user_data = sphere.user_data;
                	p.velocity += sphere.thermal_velocity(p.position);

                    add_particle_impl(data, p);
                }
            }
        }
        return ids;
	}
}
