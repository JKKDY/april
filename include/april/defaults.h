#pragma once

#include "april/monitors/monitor.h"
#include "april/monitors/output.h"
#include "april/monitors/status.h"
#include "april/monitors/performance.h"

namespace april {
	using DefaultMonitors = io::MonitorPack<io::BinaryOutput, io::ProgressBar, io::Benchmark>;

}