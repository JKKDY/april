# top level 

cmake_minimum_required(VERSION 3.28)
project(April LANGUAGES CXX)

option(APRIL_ENABLE_COVERAGE "Enable code coverage flags (disables sanitizers in Debug)" OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()


message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CMAKE_CXX_COMPILER})")



# Define the Header-Only Library
file(GLOB_RECURSE APRIL_HEADERS CONFIGURE_DEPENDS
        "${PROJECT_SOURCE_DIR}/include/april/*.h"
        "${PROJECT_SOURCE_DIR}/include/april/*.hpp"
)
add_library(April INTERFACE)  # define april as an interface library
target_sources(April INTERFACE ${APRIL_HEADERS})  # add headers to the target so IDEs can see them

# Tell targets that link to April where to find the headers
target_include_directories(April INTERFACE
        ${PROJECT_SOURCE_DIR}/include
)
target_compile_features(April INTERFACE cxx_std_23) # targets must use C++23



set_target_properties(April PROPERTIES
        CXX_STANDARD 23
        CXX_STANDARD_REQUIRED ON
)

get_target_property(APRIL_CXX_STANDARD April CXX_STANDARD)
if(APRIL_CXX_STANDARD)
    message(STATUS "C++ Standard: ${APRIL_CXX_STANDARD}")
else()
    message(STATUS "C++ Standard: Not set on target (check target_compile_features)")
endif()


# group headers in IDE
source_group(TREE "${PROJECT_SOURCE_DIR}/include" PREFIX "Header Files" FILES ${APRIL_HEADERS})


# === Warnings, optimization ===
function(set_target_build_flags target_name)
    target_compile_options(${target_name}
            PRIVATE
            # common warning flags
            $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Werror -Wpedantic>
            $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Werror -Wpedantic>
            $<$<CXX_COMPILER_ID:MSVC>:/W2 /permissive- /WX>

            # Release → -O3
            # GCC/Clang
            $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang>>:
            -O3
            -march=native           # Uses all CPU features available locally.
            -ffast-math
#            -flto                   # Link-time optimization
            -fomit-frame-pointer
            -finline-functions
            -funroll-loops
            -fno-plt
            >

            # MSVC
            $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:
            /Ox           # optimize for speed
            /GL           # whole program optimization (LTCG)
            /fp:fast      # fast-math (matches -ffast-math spirit)
            /arch:AVX2    # enable AVX2 on x64
            /favor:INTEL64
            /Oi           # intrinsic functions
            /Ot           # favor fast code
            /Ob3          # aggressive inlining
            >

            # Debug → -O0 -g
            $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU>>:-O0 -g>
            $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:Clang>>:-O0 -g>
            $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:MSVC>>:/Od /Zi>

            # Profile → -O2 -g -pg
            $<$<AND:$<CONFIG:Profile>,$<CXX_COMPILER_ID:GNU>>:-O2 -g -pg>
            $<$<AND:$<CONFIG:Profile>,$<CXX_COMPILER_ID:Clang>>:-O2 -g -pg>
    )

    # tell the linker to generate gmon.out
    target_link_options(${target_name} PRIVATE
        $<$<AND:$<CONFIG:Profile>,$<CXX_COMPILER_ID:GNU>>:-pg>
        $<$<AND:$<CONFIG:Profile>,$<CXX_COMPILER_ID:Clang>>:-pg>
    )
endfunction()


# === Sanitizer ===
function(enable_sanitizer_if_dbg target_name)
    if(APRIL_ENABLE_COVERAGE)
        message(STATUS "Skipping sanitizers for target ${target_name} because coverage is enabled.")
        return()
    endif()

    if (NOT CMAKE_BUILD_TYPE MATCHES "Debug")
        message(STATUS "Skipping AddressSanitizer for target ${target_name} because build type is not Debug.")
        return()
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
            # AppleClang does not support LeakSanitizer
            message(STATUS "Enabling AddressSanitizer + UBSan for target ${target_name} (AppleClang, Debug build)")
            target_compile_options(${target_name}
                    PRIVATE
                    -fsanitize=address,undefined -fno-omit-frame-pointer
            )
            target_link_options(${target_name}
                    PRIVATE
                    -fsanitize=address,undefined
            )
        else()
            # Linux (GCC or standard LLVM Clang)
            message(STATUS "Enabling AddressSanitizer + UBSan + Leak for target ${target_name} (GCC/Clang, Debug build)")
            target_compile_options(${target_name}
                    PRIVATE
                    -fsanitize=address,undefined,leak -fno-omit-frame-pointer
            )
            target_link_options(${target_name}
                    PRIVATE
                    -fsanitize=address,undefined,leak
            )
        endif()
    elseif (MSVC)
        message(WARNING "Skipping AddressSanitizer on MSVC for target ${target_name} to avoid LNK2038 mismatches. Enable full sanitized build instead if desired.")
    endif()
endfunction()


# === Coverage ===
function(enable_coverage target_name)
    if(NOT APRIL_ENABLE_COVERAGE)
        return()
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "Enabling LLVM coverage flags for target ${target_name}")
        target_compile_options(${target_name} PRIVATE
                -fprofile-instr-generate
                -fcoverage-mapping
        )
        target_link_options(${target_name} PRIVATE
                -fprofile-instr-generate
        )

    elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        message(STATUS "Enabling GCC coverage flags for target ${target_name}")
        target_compile_options(${target_name} PRIVATE
                --coverage
        )
        target_link_options(${target_name} PRIVATE
                --coverage
        )

    else()
        message(WARNING "Coverage is only supported on GCC and Clang. Skipping for ${target_name}.")
    endif()
endfunction()


enable_testing()

# Include subdirectories
add_subdirectory(examples)
add_subdirectory(test)
add_subdirectory(benchmark)




