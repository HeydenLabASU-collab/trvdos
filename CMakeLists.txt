cmake_minimum_required(VERSION 3.15)

project(trvdos C)

find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

add_library(vdos SHARED
    trvdos/csrc/vdos.c
)

add_library(unwrap SHARED
    trvdos/csrc/unwrap.c
)

set_target_properties(vdos PROPERTIES
    WINDOWS_EXPORT_ALL_SYMBOLS TRUE
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

set_target_properties(unwrap PROPERTIES
    WINDOWS_EXPORT_ALL_SYMBOLS TRUE
    POSITION_INDEPENDENT_CODE ON

    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/trvdos/_lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/trvdos/_lib
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/trvdos/_lib

    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/trvdos/_lib
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/trvdos/_lib
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/trvdos/_lib
)
target_link_libraries(unwrap
    PRIVATE
        GSL::gsl
        GSL::gslcblas
        OpenMP::OpenMP_C
)

install(TARGETS unwrap
    LIBRARY DESTINATION trvdos/_lib
    RUNTIME DESTINATION trvdos/_lib
    ARCHIVE DESTINATION trvdos/_lib
)