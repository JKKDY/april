# April build configuration overview
#
# Typical usage:
#   cmake -S . -B build
#   cmake --build build -j
#   ctest --test-dir build --output-on-failure
#
# Common configuration options:
#   -DAPRIL_BUILD_TESTS=ON|OFF          Build the test suite
#   -DAPRIL_BUILD_EXAMPLES=ON|OFF       Build example programs
#   -DAPRIL_ENABLE_XSIMD=ON|OFF         Use xsimd; OFF selects std::simd
#   -DAPRIL_ENABLE_OPENMP=ON|OFF        Enable the OpenMP executor backend
#   -DAPRIL_ENABLE_FAST_MATH=ON|OFF     Enable unsafe floating-point optimizations
#   -DAPRIL_MARCH=native|x86-64-v3      Select the GCC/Clang target architecture
#   -DAPRIL_MSVC_ARCH=AVX2|AVX512       Select the MSVC target instruction set
#
# Development options:
#   -DAPRIL_USE_ASAN=ON                 Enable AddressSanitizer
#   -DAPRIL_USE_UBSAN=ON                Enable UndefinedBehaviorSanitizer
#   -DAPRIL_USE_TSAN=ON                 Enable ThreadSanitizer
#   -DAPRIL_ENABLE_COVERAGE=ON          Enable coverage instrumentation
#
# Examples:
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DAPRIL_USE_ASAN=ON -DAPRIL_USE_UBSAN=ON
#   cmake -S . -B build-fast -DCMAKE_BUILD_TYPE=Release -DAPRIL_ENABLE_FAST_MATH=ON -DAPRIL_MARCH=native
#
# April is header-only. Consumers link against the April target. Tests and
# examples additionally link against April_dev_configs to receive warnings,
# optimizations, sanitizer flags, and other development-only configuration.

# This is the top level cmake
cmake_minimum_required(VERSION 3.28)


# ------------------
# PROJECT DEFINITION
# ------------------
# Deployment Target on MacOS
if(APPLE AND NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "14.0" CACHE STRING "Minimum macOS version")
endif()

project(
        April # project name: April
        VERSION 0.1.0
        DESCRIPTION "A high-performance C++23 Engine for Particle-based Simulations"
        LANGUAGES CXX # Language: C++
)

if(PROJECT_IS_TOP_LEVEL)
    message(STATUS "Configuring April v${PROJECT_VERSION} (${PROJECT_DESCRIPTION})")
endif()

# Default build type to RelWithDebInfo if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "Profile" "Assembly")


# ---------------
# COMPILE OPTIONS
# ---------------
# Component Options
option(APRIL_BUILD_TESTS      "Build the test suite" ${PROJECT_IS_TOP_LEVEL})
option(APRIL_BUILD_EXAMPLES   "Build examples" OFF)

# Backend Options (Bridges to C++ Preprocessor)
option(APRIL_ENABLE_XSIMD  "Use xsimd for vectorization" ON) # by default we use xsimd since std::simd is not mature yet on all platforms
option(APRIL_ENABLE_OPENMP "Use OpenMP executor backend" OFF)

# Code Health Options
option(APRIL_ENABLE_COVERAGE "Enable test coverage (GCC/Clang)" OFF)
option(APRIL_USE_ASAN        "Enable Address & Leak Sanitizer" OFF)
option(APRIL_USE_UBSAN       "Enable Undefined Behavior Sanitizer" OFF)
option(APRIL_USE_TSAN        "Enable Thread Sanitizer" OFF)
option(APRIL_ENABLE_STACKTRACE "Enable stack traces in APRIL assertions" OFF)

# Hardware & Optimization
option(APRIL_ENABLE_FAST_MATH "Enable unsafe floating-point optimizations (-ffast-math / /fp:fast)" OFF)
set(APRIL_MARCH "" CACHE STRING "Target architecture for GCC/Clang (e.g., native, haswell, x86-64-v3)")
set(APRIL_MSVC_ARCH "AVX2" CACHE STRING "Target instruction set for MSVC (e.g., AVX2, AVX512)")
set(APRIL_MSVC_FAVOR "NONE" CACHE STRING "Optimize for specific MSVC architecture (NONE, INTEL64, AMD64)")
set_property(CACHE APRIL_MSVC_FAVOR PROPERTY STRINGS "NONE" "INTEL64" "AMD64")

# Mutual Exclusion Checks
if(APRIL_USE_TSAN AND (APRIL_USE_ASAN OR APRIL_USE_UBSAN))
    message(FATAL_ERROR "[APRIL] ThreadSanitizer cannot be combined with AddressSanitizer or UndefinedBehaviorSanitizer.")
endif()
if(APRIL_ENABLE_COVERAGE AND (APRIL_USE_ASAN OR APRIL_USE_TSAN OR APRIL_USE_UBSAN))
    message(FATAL_ERROR "[APRIL] Coverage and sanitizers cannot be enabled together.")
endif()


# ------------------------
# APRIL LIBRARY DEFINITION
# ------------------------
# define header only library library with name April. This is what targets can link against
add_library(April INTERFACE)
target_compile_features(April INTERFACE cxx_std_23) # enforce c++23 on any target that uses April

# Public floating-point configuration
target_compile_definitions(April INTERFACE
        APRIL_FAST_MATH_ENABLED=$<BOOL:${APRIL_ENABLE_FAST_MATH}>
        APRIL_ENABLE_STACKTRACE=$<BOOL:${APRIL_ENABLE_STACKTRACE}>
)
target_compile_options(April INTERFACE
        $<$<AND:$<BOOL:${APRIL_ENABLE_FAST_MATH}>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-ffast-math>
        $<$<AND:$<BOOL:${APRIL_ENABLE_FAST_MATH}>,$<CXX_COMPILER_ID:MSVC>>:/fp:fast>
)

# Glob headers and attach them to the target for IDE
file(GLOB_RECURSE APRIL_HEADERS CONFIGURE_DEPENDS
        "${PROJECT_SOURCE_DIR}/include/april/*.h"
        "${PROJECT_SOURCE_DIR}/include/april/*.hpp"
)
target_sources(April INTERFACE ${APRIL_HEADERS})
source_group(TREE "${PROJECT_SOURCE_DIR}/include" PREFIX "Header Files" FILES ${APRIL_HEADERS})

# set the include directory
target_include_directories(April INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> # "Use this path only during the build"
        $<INSTALL_INTERFACE:include> # "Use this path after the library is installed"
)

# link to libstdc++exp if stack trace is enabled
if(APRIL_ENABLE_STACKTRACE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND
            CMAKE_SYSTEM_NAME STREQUAL "Linux")

        find_library(
                APRIL_STDCXXEXP_LIBRARY
                NAMES stdc++exp
        )

        if(NOT APRIL_STDCXXEXP_LIBRARY)
            message(FATAL_ERROR
                    "[APRIL] APRIL_ENABLE_STACKTRACE requires libstdc++exp "
                    "on Linux with GCC/libstdc++, but the library was not found."
            )
        endif()

        target_link_libraries(April INTERFACE
                ${APRIL_STDCXXEXP_LIBRARY}
        )
    endif()
endif()

# Resolve SIMD Backend
if(APRIL_ENABLE_XSIMD)
    include(FetchContent) # get xSIMD
    FetchContent_Declare(
            xsimd
            GIT_REPOSITORY https://github.com/xtensor-stack/xsimd.git
            GIT_TAG 14.0.0
            SYSTEM
    )
    FetchContent_MakeAvailable(xsimd)

    target_link_libraries(April INTERFACE xsimd)
    target_compile_definitions(April INTERFACE APRIL_SIMD_BACKEND_XSIMD)
else()
    target_compile_definitions(April INTERFACE APRIL_SIMD_BACKEND_STD_SIMD)
endif()

# Resolve Executor Backend
if(APRIL_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED) # find OpenMP
    target_link_libraries(April INTERFACE OpenMP::OpenMP_CXX)
    target_compile_definitions(April INTERFACE APRIL_EXECUTOR_BACKEND_OMP) # switch to OMP backend
endif()


# -----------
# BUILD FLAGS
# -----------
# This target holds all warnings, optimizations, and sanitizers.
# Internal targets (tests, examples) link to this. Users consuming the library do not.
add_library(April_dev_configs INTERFACE)

# Warnings
target_compile_options(April_dev_configs INTERFACE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Werror -Wpedantic>
        $<$<CXX_COMPILER_ID:MSVC>:/W2 /permissive- /WX>
)

# Optimizations
target_compile_options(April_dev_configs INTERFACE
        # Release / RelWithDebInfo (GCC/Clang/AppleClang)
        $<$<AND:$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:
        -O3
        -fomit-frame-pointer
        $<$<BOOL:${APRIL_MARCH}>:-march=${APRIL_MARCH}>
        >

        # Non-Apple GCC/Clang optimizations
        $<$<AND:$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>,$<CXX_COMPILER_ID:GNU,Clang>,$<NOT:$<PLATFORM_ID:Darwin>>>:
        -fno-plt
        >

        # Release / RelWithDebInfo (MSVC)
        $<$<AND:$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>,$<CXX_COMPILER_ID:MSVC>>:
        /O2                 # Maximize speed (replaces /Ox)
        #       /GL                 # Whole program optimization (LTCG) # commented out because of long build times
        /Gw                 # Optimize global data
        /Ob3                # Aggressive inlining (better than /Ob2)
        /Zc:inline          # Remove unreferenced COMDATs (shrinks binary/link time)
        /cgthreads8         # Use 8 threads for compilation/codegen
        /arch:${APRIL_MSVC_ARCH} # Use specified instruction set
        $<$<NOT:$<STREQUAL:${APRIL_MSVC_FAVOR},NONE>>:/favor:${APRIL_MSVC_FAVOR}>
        /bigobj
        /MP
        >

        # Debug
        $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O0 -g>
        $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:MSVC>>:/Od /Zi /bigobj>

        # Profile (GCC/Clang/AppleClang)
        $<$<AND:$<CONFIG:Profile>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-O1 -g -pg>

        # Assembly (GCC/Clang)
        # Usage: Select "Assembly" build type, then run "make <target>.s"
        $<$<AND:$<CONFIG:Assembly>,$<CXX_COMPILER_ID:GNU,Clang>,$<NOT:$<PLATFORM_ID:Darwin>>>:
        -O3
        -g
        -fverbose-asm
        -masm=intel
        $<$<BOOL:${APRIL_MARCH}>:-march=${APRIL_MARCH}>
        -fomit-frame-pointer
        -finline-functions
        -fno-plt
        >
)

# Linker flags specifically for Profile mode
target_link_options(April_dev_configs INTERFACE
        $<$<AND:$<CONFIG:Profile>,$<CXX_COMPILER_ID:GNU,Clang,AppleClang>>:-pg>
)

# Sanitizers
if(APRIL_USE_ASAN AND NOT MSVC)
    target_compile_options(April_dev_configs INTERFACE -fsanitize=address -fno-omit-frame-pointer)
    target_link_options(April_dev_configs INTERFACE -fsanitize=address)
endif()
if(APRIL_USE_UBSAN AND NOT MSVC)
    target_compile_options(April_dev_configs INTERFACE -fsanitize=undefined)
    target_link_options(April_dev_configs INTERFACE -fsanitize=undefined)
endif()
if(APRIL_USE_TSAN AND NOT MSVC)
    target_compile_options(April_dev_configs INTERFACE -fsanitize=thread -fno-omit-frame-pointer)
    target_link_options(April_dev_configs INTERFACE -fsanitize=thread)
endif()

# Coverage
if(APRIL_ENABLE_COVERAGE AND NOT MSVC)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(April_dev_configs INTERFACE -fprofile-instr-generate -fcoverage-mapping)
        target_link_options(April_dev_configs INTERFACE -fprofile-instr-generate)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(April_dev_configs INTERFACE --coverage)
        target_link_options(April_dev_configs INTERFACE --coverage)
    endif()
endif()

# handle apple compiler quirks
if(APPLE)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_definitions(April_dev_configs INTERFACE rsize_t=size_t)
        target_compile_options(April_dev_configs INTERFACE -Wno-psabi)
    endif()

    # Silence Apple Clang's warnings about features being "too new"
    if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        target_compile_options(April_dev_configs INTERFACE -Wno-unsupported-availability-guard)
    endif()
endif()

# handle MINGW compiler quirks
if(MINGW)
    # Force static linking of the GCC/Standard library to avoid DLL entry point errors
    target_link_options(April_dev_configs INTERFACE -static-libgcc -static-libstdc++ -static)
endif()


# --------------
# SUBDIRECTORIES
# --------------
if(APRIL_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(APRIL_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()


# -------------
# BUILD SUMMARY
# -------------
if(APRIL_ENABLE_XSIMD)
    set(APRIL_SIMD_DISPLAY "xsimd")
else()
    set(APRIL_SIMD_DISPLAY "std::simd")
endif()

if(MSVC)
    set(_active_arch "${APRIL_MSVC_ARCH}")
    if(NOT APRIL_MSVC_FAVOR STREQUAL "NONE")
        string(APPEND _active_arch " (Favor: ${APRIL_MSVC_FAVOR})")
    endif()
elseif(APRIL_MARCH)
    set(_active_arch "${APRIL_MARCH}")
else()
    set(_active_arch "compiler default")
endif()

message(STATUS "===========================================================")
message(STATUS "               April Configuration Summary")
message(STATUS "===========================================================")
message(STATUS "  Version:           ${PROJECT_VERSION}")
message(STATUS "  Compiler:          ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Build Type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "  Target Arch:       ${_active_arch}")
message(STATUS "  Fast Math:         ${APRIL_ENABLE_FAST_MATH}")
message(STATUS "")
message(STATUS "  SIMD Backend:      ${APRIL_SIMD_DISPLAY}")
message(STATUS "  OpenMP:            ${APRIL_ENABLE_OPENMP}")
message(STATUS "")
message(STATUS "  Sanitizers:        ASan:${APRIL_USE_ASAN} | TSan:${APRIL_USE_TSAN} | UBSan:${APRIL_USE_UBSAN}")
message(STATUS "  Coverage:          ${APRIL_ENABLE_COVERAGE}")
message(STATUS "")
message(STATUS "  Build Tests:       ${APRIL_BUILD_TESTS}")
message(STATUS "  Build Examples:    ${APRIL_BUILD_EXAMPLES}")
message(STATUS "===========================================================")