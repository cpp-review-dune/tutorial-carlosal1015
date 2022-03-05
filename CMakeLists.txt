cmake_minimum_required(VERSION 2.8.12)
project(dune-pdelab-tutorials CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# make sure we have C++14 available
dune_require_cxx_standard(MODULE "dune-pdelab-tutorials" VERSION 14)

dune_enable_all_packages()

add_subdirectory(cmake/modules)

add_subdirectory(overview)
add_subdirectory(c++)
add_subdirectory(gridinterface)
add_subdirectory(workflow)
add_subdirectory(tutorial00)
add_subdirectory(tutorial01)
add_subdirectory(tutorial02)
add_subdirectory(tutorial03)
add_subdirectory(tutorial04)
add_subdirectory(tutorial05)
add_subdirectory(tutorial06)
add_subdirectory(tutorial07)
add_subdirectory(tutorial08)
if (dune-codegen_FOUND)
  add_subdirectory(tutorial09)
endif()

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
