cmake_minimum_required(VERSION 3.20)
project(nrtk VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE)
endif()

option(NRTK_BUILD_TESTS "Build GoogleTest unit tests" ON)
option(NRTK_BUILD_RTKLIB_APPS "Build RTKLIB-EX console apps (rnx2rtkp, convbin, str2str, rtkrcv)" ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include(cmake/rtklib.cmake)
include(cmake/deps.cmake)

add_subdirectory(src)
add_subdirectory(apps)
if(NRTK_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
