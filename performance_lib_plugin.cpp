/**
 * performance_lib plugin
 */

#include "bridge/util.h"


using namespace std;
using namespace ug::bridge;

void SetCPUFrequency(double f)
{
	UG_LOG("Setting CPU Freq to " << f << "\n");
}

/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_PerformanceLib(Registry* reg, string grp)
{
	reg->add_function("SetCPUFrequency", SetCPUFrequency, grp, "", "frequency", "set the cpu frequency to f");
}// namespace ug
