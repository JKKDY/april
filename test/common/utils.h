#pragma once

#include "april/april.hpp"

using namespace april;


template<core::IsSystem System>
System::ParticleRec get_particle(System& sys, size_t index) {
	constexpr auto all_fields = ParticleField::all;

	auto p_ref = sys.template at<all_fields>(index);

	typename System::ParticleRec rec;
	rec.id          = p_ref.id;
	rec.type        = p_ref.type;
	rec.position    = p_ref.position;
	rec.velocity    = p_ref.velocity;
	rec.force       = p_ref.force;
	rec.old_position = p_ref.old_position;
	rec.state       = p_ref.state;
	rec.mass        = p_ref.mass;
	rec.attributes   = p_ref.attributes;

	return rec;
}

template<core::IsSystem System>
System::ParticleRec get_particle_by_id(System& sys, ParticleID id) {
	constexpr auto all_fields = ParticleField::all;

	auto p_ref = sys.template view_id<all_fields>(id);

	typename System::ParticleRec rec;
	rec.id          = p_ref.id;
	rec.type        = p_ref.type;
	rec.position    = p_ref.position;
	rec.velocity    = p_ref.velocity;
	rec.force       = p_ref.force;
	rec.old_position = p_ref.old_position;
	rec.state       = p_ref.state;
	rec.mass        = p_ref.mass;
	rec.attributes   = p_ref.attributes;

	return rec;
}

template<core::IsSystem System>
std::vector<typename System::ParticleRec> export_particles(System& sys) {
	std::vector<typename System::ParticleRec> records;
	records.reserve(sys.size());

	sys.for_each_particle(april::scalar_kernel<ParticleField::none>(
		[&](size_t idx, auto &&) {
			records.push_back(get_particle(sys, idx));
		}
	));

	return records;
}



template<core::IsSystem System>
void simulate_single_step(System& sys) {
	constexpr ParticleField edit_fields = ParticleField::old_position | ParticleField::position | ParticleField::velocity;

	for (size_t pid = sys.min_id(); pid < sys.max_id(); ++pid) {
		if (!sys.contains_id(pid)) continue;
		auto p = sys.template at_id<edit_fields>(pid);
		p.old_position = p.position;
		p.position = p.old_position + p.velocity; // simulate one step
	}
}


inline Particle make_particle(
	const ParticleType type, const vec3 & position, const vec3 & velocity, const double mass,
	const ParticleState state = ParticleState::ALIVE, const std::optional<ParticleID> id = std::nullopt) {

	Particle p;
	p.type = type;
	p.position = position;
	p.velocity = velocity;
	p.mass = mass;
	p.state = state;
	p.id = id;
	return p;
}


template<ParticleField Read, ParticleField Write = Read, typename RecordT>
auto make_particle_source(RecordT& record) {
	auto get_field = [&]<ParticleField F>() {
		if constexpr (F == ParticleField::position) return &record.position;
		else if constexpr (F == ParticleField::velocity) return &record.velocity;
		else if constexpr (F == ParticleField::force) return &record.force;
		else if constexpr (F == ParticleField::old_position) return &record.old_position;
		else if constexpr (F == ParticleField::mass) return &record.mass;
		else if constexpr (F == ParticleField::state) return &record.state;
		else if constexpr (F == ParticleField::type) return &record.type;
		else if constexpr (F == ParticleField::id) return &record.id;
		else if constexpr (F == ParticleField::attributes) return &record.attributes;
	};

	return particle::internal::make_particle_source<Read, Write>(get_field);
}
