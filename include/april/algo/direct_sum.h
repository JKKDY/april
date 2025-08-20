#pragma once


#include "april/algo/container.h"

namespace april::core::impl {
	class DirectSum;
}

namespace april::core {

	class DirectSum {
		using Container = impl::DirectSum;
	};

	namespace impl {
		class DirectSum final : public Container {
		public:
			using Container::Container;

			void build() override;
			void calculate_forces() override;
		};
	}
}