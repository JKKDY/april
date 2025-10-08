#pragma once

#include "april/monitors/monitor.h"
#include "april/monitors/binary_output.h"
#include "april/monitors/progressbar.h"
#include "april/monitors/benchmark.h"

namespace april {
	using DefaultMonitors = monitor::MonitorPack<monitor::BinaryOutput, monitor::ProgressBar, monitor::Benchmark>;

}