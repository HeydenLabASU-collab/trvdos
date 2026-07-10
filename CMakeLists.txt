cmake_minimum_required(VERSION 3.15)

project(trvdos C)

message(STATUS "CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
message(STATUS "CMAKE_SYSTEM_PREFIX_PATH=${CMAKE_SYSTEM_PREFIX_PATH}")

find_package(PkgConfig)

if(PkgConfig_FOUND)
    pkg_check_modules(GSL gsl)

    message(STATUS "PkgConfig GSL_FOUND=${GSL_FOUND}")
    message(STATUS "GSL_INCLUDE_DIRS=${GSL_INCLUDE_DIRS}")
    message(STATUS "GSL_LIBRARY_DIRS=${GSL_LIBRARY_DIRS}")
    message(STATUS "GSL_LIBRARIES=${GSL_LIBRARIES}")
endif()

find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

add_library(vdos SHARED
    trvdos/csrc/vdos.c
)

target_link_libraries(vdos PRIVATE
    GSL::gsl
    GSL::gslcblas
    OpenMP::OpenMP_C
)

set_target_properties(vdos PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/trvdos/_lib
)

install(TARGETS vdos
    LIBRARY DESTINATION trvdos/_lib
    RUNTIME DESTINATION trvdos/_lib
    ARCHIVE DESTINATION trvdos/_lib
)