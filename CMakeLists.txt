################################################################################
# PerformanceLib
#
# The PerformanceLib plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

# Name of your plugin, path to your plugin (relative to ug's root path) and sources.
set(pluginName	PerformanceLib)
set(pluginPath	plugins/experimental/template)


project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")




#### here go all sources
set(SOURCES		performance_lib_plugin.cpp)



################################################################################
# The code below doesn't have to be changed (usually)
################################################################################


if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${pluginPath} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
endif(buildEmbeddedPlugins)

