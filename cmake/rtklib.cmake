# RTKLIB-EX (rtklibexplorer) built from the pinned git submodule.
#
# We deliberately do not add_subdirectory() the upstream CMakeLists.txt:
#  * it declares NFREQ/NEXOBS/ENA* with add_definitions(), i.e. only for its own
#    directory. Code that includes rtklib.h from outside would silently see a
#    different obsd_t/rtcm_t layout. Here the definitions are PUBLIC on the
#    target, so every consumer compiles against the same struct layout.
#  * it writes libraries and executables into the submodule source tree.
# The upstream sources are used unmodified.

set(RTKLIB_DIR ${PROJECT_SOURCE_DIR}/third_party/RTKLIB)
if(NOT EXISTS ${RTKLIB_DIR}/src/rtklib.h)
  message(FATAL_ERROR "RTKLIB-EX submodule not found. Run: git submodule update --init")
endif()

file(STRINGS ${RTKLIB_DIR}/src/rtklib.h _rtklib_patch REGEX "^#define PATCH_LEVEL")
string(REGEX MATCH "\"([^\"]*)\"" _ "${_rtklib_patch}")
set(RTKLIB_VERSION "${CMAKE_MATCH_1}")
message(STATUS "RTKLIB-EX version: ${RTKLIB_VERSION}")

# same source set as upstream src/CMakeLists.txt (comnav.c/tersus.c are not
# maintained upstream and do not compile against the current headers)
file(GLOB RTKLIB_SRCS ${RTKLIB_DIR}/src/*.c)
foreach(rcv binex crescent javad novatel nvs rt17 septentrio skytraq swiftnav ublox unicore)
  list(APPEND RTKLIB_SRCS ${RTKLIB_DIR}/src/rcv/${rcv}.c)
endforeach()

add_library(rtklib_ex STATIC ${RTKLIB_SRCS})
target_include_directories(rtklib_ex SYSTEM PUBLIC ${RTKLIB_DIR}/src)
target_compile_definitions(rtklib_ex PUBLIC
  ENAGLO ENAQZS ENACMP ENAGAL ENAIRN NFREQ=3 NEXOBS=3 TRACE SVR_REUSEADDR)
set_target_properties(rtklib_ex PROPERTIES C_STANDARD 99 POSITION_INDEPENDENT_CODE ON)
target_compile_options(rtklib_ex PRIVATE -w)
find_package(Threads REQUIRED)
target_link_libraries(rtklib_ex PUBLIC Threads::Threads m)

if(NRTK_BUILD_RTKLIB_APPS)
  set(_consapp ${RTKLIB_DIR}/app/consapp)
  foreach(app rnx2rtkp convbin str2str pos2kml)
    add_executable(${app} ${_consapp}/${app}/${app}.c)
    target_link_libraries(${app} PRIVATE rtklib_ex)
    target_compile_options(${app} PRIVATE -w)
  endforeach()
  if(UNIX)
    add_executable(rtkrcv ${_consapp}/rtkrcv/rtkrcv.c ${_consapp}/rtkrcv/vt.c)
    target_include_directories(rtkrcv PRIVATE ${_consapp}/rtkrcv)
    target_link_libraries(rtkrcv PRIVATE rtklib_ex)
    target_compile_options(rtkrcv PRIVATE -w)
  endif()
endif()
