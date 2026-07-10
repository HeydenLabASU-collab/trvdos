cmake_minimum_required(VERSION 3.15)

project(trvdos C)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GSL REQUIRED gsl)

find_package(OpenMP REQUIRED)

add_library(vdos SHARED
    trvdos/csrc/vdos.c
)

target_include_directories(vdos PRIVATE
    ${GSL_INCLUDE_DIRS}
)

target_link_directories(vdos PRIVATE
    ${GSL_LIBRARY_DIRS}
)

target_link_libraries(vdos PRIVATE
    ${GSL_LIBRARIES}
    OpenMP::OpenMP_C
)

target_compile_options(vdos PRIVATE
    ${GSL_CFLAGS_OTHER}
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