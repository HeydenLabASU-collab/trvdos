cmake_minimum_required(VERSION 3.15)

project(mda_3d_2pt C)

find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

add_library(vdos SHARED
    mda_3d_2pt/csrc/vdos.c
)

set_target_properties(vdos PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

target_link_libraries(vdos
    PRIVATE
    GSL::gsl
    GSL::gslcblas
    OpenMP::OpenMP_C)

set_target_properties(vdos PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY
        ${CMAKE_SOURCE_DIR}/mda_3d_2pt/_lib
)
