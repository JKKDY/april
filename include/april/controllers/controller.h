#pragma once

namespace april::env {

	class Controller {
	public:
		explicit Controller(const size_t call_frequency) : call_frequency_m(call_frequency) {}

		[[nodiscard]] size_t call_frequency() const { return call_frequency_m; }

		void apply(this auto&& self, size_t step) {
			// static_assert(
			// 	requires { self.write_output(step); },
			// 	"OutputWriter requires a write_output(size_t, const std::vector<Particle>&) method"
			// );
			// self.write_output(step);
		}

	private:
		size_t call_frequency_m{};

	};
}