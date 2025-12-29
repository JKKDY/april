#pragma once

# include "april/april.h"

using namespace april;


template<core::IsSystem System>
typename System::ParticleRec get_particle(System& sys, size_t index) {
	constexpr auto all_fields = to_field_mask(env::Field::all);

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
	rec.user_data   = p_ref.user_data;

	return rec;
}

template<core::IsSystem System>
typename System::ParticleRec get_particle_by_id(System& sys, ParticleID id) {
	constexpr auto all_fields = to_field_mask(env::Field::all);

	auto p_ref = sys.template at_id<all_fields>(id);

	typename System::ParticleRec rec;
	rec.id          = p_ref.id;
	rec.type        = p_ref.type;
	rec.position    = p_ref.position;
	rec.velocity    = p_ref.velocity;
	rec.force       = p_ref.force;
	rec.old_position = p_ref.old_position;
	rec.state       = p_ref.state;
	rec.mass        = p_ref.mass;
	rec.user_data   = p_ref.user_data;

	return rec;
}

template<core::IsSystem System>
std::vector<typename System::ParticleRec> export_particles(System& sys) {
	std::vector<typename System::ParticleRec> records;
	records.reserve(sys.size());

	for (size_t i = 0; i < sys.size(); ++i) {
		records.push_back(get_particle(sys, i));
	}

	return records;
}



template<core::IsSystem System>
void simulate_single_step(System& sys) {
	constexpr env::FieldMask edit_fields = env::Field::old_position | env::Field::position | env::Field::velocity;

	for (size_t pid = 0; pid < sys.size(); ++pid) {
		auto p = sys.template at<edit_fields>(pid);
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
