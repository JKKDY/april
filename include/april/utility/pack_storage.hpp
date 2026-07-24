#pragma once
#include <tuple>
#include <vector>

#include "april/base/concepts.hpp"

namespace april::utility::internal {


	template<typename... Ts>
	inline constexpr bool unique_types_v = true;

	template<typename T, typename... Ts>
	inline constexpr bool unique_types_v<T, Ts...> =
		(!same_as_any<T, Ts...>) &&
		unique_types_v<Ts...>;



	template<class... Ts>
	class PackStorage {
		static_assert(
			unique_types_v<Ts...>,
			"PackStorage component types must be unique."
		);
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

		template<typename Func>
		void for_each_list(Func&& function) const {
			std::apply(
				[&](const auto&... lists) {
					(function(lists), ...);
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

		template<typename Func>
		void for_each_item(Func&& f) const {
			for_each_list([&](const auto& list) {
				for (const auto& item : list) {
					f(item);
				}
			});
		}

		// Get the vector for a specific type T
		template<typename T>
		std::vector<T>& get_list() {
			return std::get<std::vector<T>>(components);
		}

		// remove all components
		void clear() {
			std::apply(
			   [](auto&... list) {
				  (list.clear(), ...);
			   },
			   components
			);
		}

	private:
		std::tuple<std::vector<Ts>...> components{};
	};
}















