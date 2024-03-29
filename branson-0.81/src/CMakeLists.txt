################################################################################
# NOTES ON COMPILATION
################################################################################

# Load the modules for each libary and the environment variables will be set
# correctly (load cmake, compiler, mpi, metis, parmetis and optionally
# hdf5 and silo)

# Enbable debug mode by passing -DCMAKE_BUILD_TYPE=Debug to CMake, default is
# Release

cmake_minimum_required (VERSION 3.11)
project (BRANSON
  VERSION 0.8
  DESCRIPTION "Branson can be used to study different algorithms for parallel Monte Carlo transport. Currently it contains particle passing and mesh passing methods for domain decomposition."
#  HOMEPAGE URL "https://github.com/lanl/branson" # needs cmake 3.12+
  LANGUAGES CXX C)

# Build system support files are located here.
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config )

message( "
This is ${CMAKE_PROJECT_NAME} v. ${PROJECT_VERSION}.
")

#------------------------------------------------------------------------------#
# User configurable options

set (BRANSON_N_GROUPS ${N_GROUPS})
if (NOT BRANSON_N_GROUPS)
  message(STATUS
    "HARDCODED NUMBER OF GROUPS NOT SET, DEFAULTING TO ONE GROUP (GRAY)
     use '-DN_GROUPS=<int>' when running cmake to select multigroup.\n")
  set(BRANSON_N_GROUPS "1")
else ()
  message(STATUS
    "HARDCODED NUMBER OF GROUPS SET TO ${BRANSON_N_GROUPS} (FAUX MULTIGROUP)\n")
endif()
set(BRANSON_N_GROUPS ${BRANSON_N_GROUPS} CACHE STRING "Number of energy groups." FORCE)

#------------------------------------------------------------------------------#
# Record the site name

site_name( SITENAME )
string( REGEX REPLACE "([A-z0-9]+).*" "\\1" SITENAME ${SITENAME} )
include (lanl-setup)
set( SITENAME ${SITENAME} CACHE "STRING" "Name of the current machine" FORCE)

#------------------------------------------------------------------------------#
# Setup compiler options

# C11 support:
set( CMAKE_C_STANDARD 11 )

# C++14 support:
set( CMAKE_CXX_STANDARD 11 )
set( CXX_STANDARD_REQUIRED ON )

# Do not enable extensions (e.g.: --std=gnu++11)
set( CMAKE_CXX_EXTENSIONS OFF )
set( CMAKE_C_EXTENSIONS   OFF )

#------------------------------------------------------------------------------#
# Build type and custom compiler flags

if ( "${CMAKE_BUILD_TYPE}notset" STREQUAL "notset" AND
    NOT DEFINED CMAKE_CONFIGURATION_TYPES )
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Release|Debug|RelWithDebInfo"
    FORCE)
endif ()

# GCC options for address or undefined sanitizing
#set(GCC_SANITIZE_COMPILE_FLAGS "-fsanitize=address")
#set(GCC_SANITIZE_LINK_FLAGS "-fsanitize=address")
#set(GCC_SANITIZE_COMPILE_FLAGS "-fsanitize=undefined")
#set(GCC_SANITIZE_LINK_FLAGS "-fsanitize=undefined")

# add compiler flags
if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
  set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0")
  # optionally add sanitize flags
  #set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${GCC_SANITIZE_COMPILE_FLAGS}")
  #set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_SANITIZE_LINK_FLAGS}")
elseif( CMAKE_CXX_COMPILER_ID STREQUAL "Intel" )
  set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -O0")
elseif( CMAKE_CXX_COMPILER_ID STREQUAL "Cray" )
  set( CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -DR123_USE_GNU_UINT128=0")
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DR123_USE_GNU_UINT128=0")
endif()

#------------------------------------------------------------------------------#
# Look for Third Party Libraries (metis, parmetis, etc.)
add_subdirectory(pugixml)

include(find_tpls)
setupTPLs()
#------------------------------------------------------------------------------#
# Report build configuration

# Summary of TPLs
include(FeatureSummary)
message(" ")
feature_summary( WHAT ALL
  INCLUDE_QUIET_PACKAGES
  FATAL_ON_MISSING_REQUIRED_PACKAGES
  QUIET_ON_EMPTY )

# Summary of user-selectable build options
message( "\nBuild Summary:\n")
message( STATUS "Machine name    : ${SITENAME}")
report_lanl_hpc_features()
if( ${BRANSON_N_GROUPS} GREATER 1 )
  message( STATUS "Energy groups   : ${BRANSON_N_GROUPS}")
else()
  message( STATUS "Energy groups   : Gray (1-group)")
endif()
message( STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message( STATUS "Compiler        : ${CMAKE_CXX_COMPILER}")
if( CMAKE_CONFIGURATION_TYPES )
  message( STATUS "Compiler Flags (All)    : ${CMAKE_CXX_FLAGS}")
  message( STATUS "Compiler Flags (Debug)  : ${CMAKE_CXX_FLAGS_DEBUG}")
  message( STATUS "Compiler Flags (Release): ${CMAKE_CXX_FLAGS_RELEASE}")
else()
  string( TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER )
  message( STATUS "Compiler Flags  : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")
endif()
message("\n")

#------------------------------------------------------------------------------#
# Set up libraries for Cray RMA routines

# set(DMAPP_DYNAMIC -Wl,--whole-archive,-ldmapp,--no-whole-archive)

#------------------------------------------------------------------------------#
# Generate config.h

configure_file(config.h.in ${PROJECT_BINARY_DIR}/config.h)

#------------------------------------------------------------------------------#
# Targets

add_executable(BRANSON main.cc)
target_include_directories( BRANSON PRIVATE
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> ${PROJECT_SOURCE_DIR}/pugixml/src/)
set( branson_deps
  ParMETIS::parmetis
  MPI::MPI_CXX
  MPI::MPI_C
  METIS::metis 
  pugixml)
if( VIZ_LIBRARIES_FOUND )
  set( branson_deps "Silo::silo;${HDF5_LIBRARIES};${branson_deps}" )
endif()
target_link_libraries( BRANSON PUBLIC ${branson_deps} )

#------------------------------------------------------------------------------#
# Testing

enable_testing()
add_subdirectory(test)

#------------------------------------------------------------------------------#
# Targets for installation

install(TARGETS BRANSON DESTINATION bin)

#------------------------------------------------------------------------------#
# End src/CMakeLists.txt
#------------------------------------------------------------------------------#
