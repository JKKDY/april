#pragma once
#include <tuple>
#include <vector>

#include "april/common.h"


namespace april::shared::internal {

	template<typename P>
	concept IsPack = requires {
		typename P::is_pack_tag;
		typename std::tuple_element_t<0, typename P::types>;
	};

	template<class... Ts>
	class PackStorage {
	public:
		PackStorage() = default;

		// Add one component of type T
		template <class T>
		requires (same_as_any<T, Ts...>)
		void add(T component) {
			std::get<std::vector<T>>(components).push_back(std::move(component));
		}

		// Add multiple at once
		template <class... Cs>
		void add_many(Cs&&... comps) {
			(add(std::forward<Cs>(comps)), ...);
		}

		// Generic loop over all component vectors
		// e.g., for_each_list([&](auto& vec){ for(auto& item : vec) { ... } });
		template<typename Func>
		void for_each_list(Func&& f) {
			std::apply(
				[&](auto&... list) {
					(f(list), ...);
				},
				components
			);
		}

		// Generic loop over every single component item
		// e.g., for_each_item([&](auto& item){ item.update(...); });
		template<typename Func>
		void for_each_item(Func&& f) {
			for_each_list([&](auto& list) {
				for (auto& item : list) {
					f(item);
				}
			});
		}

		// Get the vector for a specific type T
		template<typename T>
		std::vector<T>& get_list() {
			return std::get<std::vector<T>>(components);
		}

	private:
		std::tuple<std::vector<Ts>...> components{};
	};
}
