# top level 

cmake_minimum_required(VERSION 3.15)
project(April LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# For all targets in this directory and below
if (MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(GLOBAL_OPT_FLAGS -O3)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(GLOBAL_OPT_FLAGS -O0 -g)
endif()

add_compile_options(${GLOBAL_OPT_FLAGS})

# Include subdirectories
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(test)


