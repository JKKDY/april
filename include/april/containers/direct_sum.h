#pragma once


#include "april/containers/container.h"

namespace april::core {

	class DirectSum final : public Container {
	public:
		using Container::Container;

		void build() override;
		void calculate_forces() override;
	};
}