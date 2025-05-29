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

# Include subdirectories
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(test)


