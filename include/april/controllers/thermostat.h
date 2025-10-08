#pragma once

#include <limits>

namespace april::controller {
    class Environment;

    struct Thermostat {
        double init_temp;           ///< The initial temperature.
        double target_temp;         ///< The target temperature.
        double max_temp_change;     ///< The maximal absolute temperature change allowed for one application of the thermostat.
		unsigned int call_frequency;///< The frequency of thermostat calls.
    };

    namespace impl
    {
      
        class Thermostat {
        public:
            static constexpr double INF_TEMP = std::numeric_limits<double>::max();
            static constexpr double DONT_CARE_TEMP = -1;

            explicit Thermostat(double init_T = DONT_CARE_TEMP, double target_T = DONT_CARE_TEMP, double dT = INF_TEMP);
            void init(double init_T, double target_T, double dT);
            void set_initial_temperature(Environment& env) const;
            void adjust_temperature(Environment& env) const;

        private:
            double init_temp;        ///< The initial temperature.
            double target_temp;      ///< The target temperature.
            double max_temp_change;  ///< The maximal absolute temperature change allowed for one application of the thermostat.
        };
    }
} // namespace md::env



