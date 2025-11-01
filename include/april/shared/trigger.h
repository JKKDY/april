#pragma once

#include <cstddef>
#include <concepts>

#include "april/core/context.h"

namespace april::trigger {

	struct Trigger {
		template<class S>
		bool operator()(this const auto&& self, const core::SystemContext<S>& sys) const noexcept {
			static_assert(
				requires { {self.should_fire(sys) } -> std::same_as<bool>; },
				"Trigger subclass must implement operator()"
			);
			return self.should_fire(sys);
		}
	};

	template <class T>
	concept IsTrigger = std::derived_from<T, Trigger>;



	// ---- Pre defined Triggers ----

	// -- step based triggers: --
	// trigger every N steps
	struct Every : Trigger {
		explicit Every(const size_t N, const size_t offset = 0):
			N(N), offset(offset) {}

		template<class S>
		bool should_fire(const core::SystemContext<S>& sys) {
			return (sys.step() + offset) % N == 0;
		}

	private:
		size_t N, offset;
	};

	// start triggering after *step*
	struct After : Trigger {
		explicit After(const size_t step): step(step) {}

		template<class S>
		bool should_fire(const core::SystemContext<S>& sys) {
			return (sys.step() >= step);
		}

	private:
		size_t step;
	};

	// Trigger while current step ∈ [start, end)
	struct Between : Trigger {
		Between(const size_t start, const size_t end) : start(start), end(end) {}

		template<class S>
		bool should_fire(const core::SystemContext<S>& sys) const noexcept {
			return sys.step() >= start && sys.step() < end;
		}

	private:
		size_t start, end;
	};

	// Trigger exactly at one step
	struct AtStep : Between {
		explicit AtStep(const size_t step) : Between(step, step+1) {}
	};


	// -- Time-based triggers --
	// Trigger periodically every `period` time units, with optional offset
	struct Periodically : Trigger {
		explicit Periodically(const double period, const double offset = 0.0)
			: period(period), offset(offset), last(offset - period) {}

		template<class S>
		bool should_fire(const core::SystemContext<S>& sys) noexcept {
			if (sys.time() - last >= period) {
				last = sys.time();
				return true;
			}
			return false;
		}

	private:
		double period, offset, last;
	};

	// Trigger only after a given time
	struct AfterTime : Trigger {
		explicit AfterTime(const double t) : t(t) {}

		template<class S>
		bool should_fire(const core::SystemContext<S>& sys) const noexcept {
			return sys.time() >= t;
		}

	private:
		double t;
	};

	// Trigger only while time ∈ [t_start, t_end)
	struct BetweenTime : Trigger {
		BetweenTime(const double t_start, const double t_end)
			: t_start(t_start), t_end(t_end) {}

		template<class S>
		bool should_fire(const core::SystemContext<S>& sys) const noexcept {
			return sys.time() >= t_start && sys.time() < t_end;
		}

	private:
		double t_start, t_end;
	};


	// -- Generic triggers --
	// Always true
	struct Always : Trigger {
		template<class S>
		bool should_fire(const core::SystemContext<S>&) const noexcept {
			return true;
		}
	};



	// ---- Logical Operators ----
	template<class L, class R>
	struct AndTrigger : Trigger<AndTrigger<L,R>> {
		L lhs;
		R rhs;

		template<class S>
		bool should_fire(const core::SystemContext<S>& sys) const noexcept {
			return lhs(sys) && rhs(sys);
		}
	};

	template<class L, class R>
	struct OrTrigger : Trigger<OrTrigger<L,R>> {
		L lhs;
		R rhs;

		template<class S>
		bool should_fire(const core::SystemContext<S>& sys) const noexcept {
			return lhs(sys) || rhs(sys);
		}
	};

	template<class T>
	struct NotTrigger : Trigger<NotTrigger<T>> {
		T inner;

		template<class S>
		bool should_fire(const core::SystemContext<S>& sys) const noexcept {
			return !inner(sys);
		}
	};

	template<class L, class R>
	constexpr auto operator&&(const Trigger<L>& lhs, const Trigger<R>& rhs) noexcept {
		return AndTrigger<L, R>{static_cast<const L&>(lhs), static_cast<const R&>(rhs)};
	}

	template<class L, class R>
	constexpr auto operator||(const Trigger<L>& lhs, const Trigger<R>& rhs) noexcept {
		return OrTrigger<L, R>{static_cast<const L&>(lhs), static_cast<const R&>(rhs)};
	}

	template<class T>
	constexpr auto operator!(const Trigger<T>& t) noexcept {
		return NotTrigger<T>{static_cast<const T&>(t)};
	}
}
