#pragma once

namespace april::controller::impl {

	class Controller {
	public:
		explicit Controller(const size_t call_frequency) : call_frequency_m(call_frequency) {}

		[[nodiscard]] size_t call_frequency() const { return call_frequency_m; }

		void apply(this auto&& self, size_t step) {

		}

	private:
		size_t call_frequency_m{};
	};
}