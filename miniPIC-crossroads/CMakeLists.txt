cmake_minimum_required (VERSION 2.8.11)

SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})

FIND_PACKAGE(Trilinos REQUIRED)
project (mini-PIC CXX)

# Make sure to use same compilers and flags as Trilinos
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
INCLUDE_DIRECTORIES (${Trilinos_INCLUDE_DIRS} )

option (ENABLE_PIC_TESTS "enable tests" OFF)

add_subdirectory(src)

IF(ENABLE_PIC_TESTS) 
enable_testing()
add_subdirectory(test)
ENDIF()

