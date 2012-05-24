################################################################################
# TEMPLATE
#
# The template plugin.
################################################################################

cmake_minimum_required(VERSION 2.6)

project(TEMPLATE)

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")


# sources
set(SOURCES	template_plugin.cpp)



if(STATIC)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/template ${SOURCES})
else(STATIC)
	add_library(Template SHARED ${SOURCES})
	target_link_libraries (Template ug4)
endif(STATIC)
