cmake_minimum_required(VERSION 3.15)

project(trvdos C)

find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

add_library(vdos SHARED
    trvdos/csrc/vdos.c
)

set_target_properties(vdos PROPERTIES
    POSITION_INDEPENDENT_CODE ON

    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/trvdos/_lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/trvdos/_lib
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/trvdos/_lib

    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/trvdos/_lib
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/trvdos/_lib
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/trvdos/_lib
)
target_link_libraries(vdos
    PRIVATE
        GSL::gsl
        GSL::gslcblas
        OpenMP::OpenMP_C
)

install(TARGETS vdos
    LIBRARY DESTINATION trvdos/_lib
    RUNTIME DESTINATION trvdos/_lib
    ARCHIVE DESTINATION trvdos/_lib
)