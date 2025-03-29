# top level 

cmake_minimum_required(VERSION 3.15)
project(April LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include subdirectories
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(test)