# Third-party C++ dependencies, pinned. A system package is used when present,
# otherwise the source is fetched into the build tree (nothing is installed).
include(FetchContent)

find_package(yaml-cpp 0.8 QUIET CONFIG)
if(NOT yaml-cpp_FOUND)
  set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
  set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
  set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
  set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG        yaml-cpp-0.9.0
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(yaml-cpp)
endif()
if(NOT TARGET yaml-cpp::yaml-cpp)
  add_library(yaml-cpp::yaml-cpp ALIAS yaml-cpp)
endif()

if(NRTK_BUILD_TESTS)
  find_package(GTest 1.14 QUIET CONFIG)
  if(NOT GTest_FOUND)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(googletest
      GIT_REPOSITORY https://github.com/google/googletest.git
      GIT_TAG        v1.17.0
      GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(googletest)
  endif()
endif()
