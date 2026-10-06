file(GLOB NRTK_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
if(NRTK_TEST_SOURCES)
  add_executable(nrtk_tests ${NRTK_TEST_SOURCES})
  target_link_libraries(nrtk_tests PRIVATE nrtk GTest::gtest_main)
  target_compile_definitions(nrtk_tests PRIVATE NRTK_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  include(GoogleTest)
  gtest_discover_tests(nrtk_tests DISCOVERY_TIMEOUT 60)
endif()

# end-to-end regression (nrtk-sim -> nrtk-gen -> obscmp/rnx2rtkp)
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND AND NRTK_BUILD_RTKLIB_APPS)
  add_test(NAME regression_sim_loo
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/regression/run_regression.py
            ${CMAKE_BINARY_DIR}/bin ${PROJECT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/regression)
  set_tests_properties(regression_sim_loo PROPERTIES TIMEOUT 900 LABELS regression)
endif()
