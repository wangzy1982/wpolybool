cmake_minimum_required(VERSION 3.15)
project(wpolybool_project LANGUAGES CXX C)

add_library(wpolybool SHARED)

target_include_directories(wpolybool
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdpart/wsolver/include
)

file(GLOB_RECURSE WPOLY_SRC_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c
)
target_sources(wpolybool PRIVATE ${WPOLY_SRC_FILES})

if(WIN32)
    set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "Supported build configurations" FORCE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)

    target_compile_definitions(wpolybool
        PRIVATE
            _WINDOWS
            WPOLY_EXPORTS
    )
    
    target_link_directories(wpolybool
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/thirdpart/wsolver/lib
    )
    target_link_libraries(wpolybool
        PRIVATE
            wsolver.lib
    )
endif()

set_target_properties(wpolybool PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib
)
set_target_properties(wpolybool PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/bin
    ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/lib
)