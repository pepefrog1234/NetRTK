# libnrtk: network RTK core library (C11 / C++17)
set(NRTK_C_SOURCES
  common/common.c
  common/tropo.c
  common/gmf.c
  common/antenna.c
  common/prep.c
)
set(NRTK_CXX_SOURCES)

foreach(f IN ITEMS
    ingest/ingest.c
    netpos/netpos.c
    vrsgen/vrsgen.c
    rtcm_out/rtcm_out.c
    datum/datum.c)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${f})
    list(APPEND NRTK_C_SOURCES ${f})
  endif()
endforeach()
foreach(f IN ITEMS
    cellmgr/cellmgr.cpp
    cellmgr/delaunay.cpp
    atmos/atmos.cpp
    epoch_sync/epoch_sync.cpp
    config/config.cpp
    monitor/monitor.cpp
    ntrip_front/ntrip_server.cpp
    ntrip_front/gga.cpp
    pipeline/network.cpp
    pipeline/offline.cpp
    pipeline/realtime.cpp)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${f})
    list(APPEND NRTK_CXX_SOURCES ${f})
  endif()
endforeach()

add_library(nrtk STATIC ${NRTK_C_SOURCES} ${NRTK_CXX_SOURCES})
target_include_directories(nrtk PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(nrtk PUBLIC rtklib_ex yaml-cpp::yaml-cpp)
target_compile_options(nrtk PRIVATE
  $<$<COMPILE_LANGUAGE:C>:-Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers>
  $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers>)
