#pragma once

#include "april/io/monitor.h"
#include "april/io/output.h"
#include "april/io/status.h"
#include "april/io/performance.h"

namespace april {
	using DefaultMonitors = io::MonitorPack<io::BinaryOutput, io::ProgressBar, io::Benchmark>;

}