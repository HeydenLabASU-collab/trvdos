if(DEFINED ENV{CONDA_PREFIX})
    list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif()

cmake_minimum_required(VERSION 3.15)

project(trvdos C)

message(STATUS "CONDA_PREFIX=$ENV{CONDA_PREFIX}")
message(STATUS "GSL_ROOT_DIR=${GSL_ROOT_DIR}")
find_program(GSL_CONFIG_EXECUTABLE gsl-config)
message(STATUS "gsl-config=${GSL_CONFIG_EXECUTABLE}")

find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

add_library(vdos SHARED
    trvdos/csrc/vdos.c)

set_target_properties(vdos PROPERTIES
    POSITION_INDEPENDENT_CODE ON)

target_link_libraries(vdos
    PRIVATE
    GSL::gsl
    GSL::gslcblas
    OpenMP::OpenMP_C)

set_target_properties(vdos PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY
        ${CMAKE_SOURCE_DIR}/trvdos/_lib)

install(TARGETS vdos
        LIBRARY DESTINATION trvdos/_lib
        RUNTIME DESTINATION trvdos/_lib
        ARCHIVE DESTINATION trvdos/_lib)