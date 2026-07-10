cmake_minimum_required(VERSION 3.15)

project(trvdos C)

find_program(GSL_CONFIG REQUIRED gsl-config)

execute_process(
    COMMAND ${GSL_CONFIG} --prefix
    OUTPUT_VARIABLE GSL_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "GSL prefix = ${GSL_PREFIX}")

include_directories(${GSL_PREFIX}/include)
link_directories(${GSL_PREFIX}/lib)

find_package(OpenMP REQUIRED)

add_library(vdos SHARED
    trvdos/csrc/vdos.c
)

set_target_properties(vdos PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/trvdos/_lib
)

target_link_libraries(vdos
    PRIVATE
        gsl
        gslcblas
        OpenMP::OpenMP_C
)

install(TARGETS vdos
    LIBRARY DESTINATION trvdos/_lib
    RUNTIME DESTINATION trvdos/_lib
    ARCHIVE DESTINATION trvdos/_lib
)