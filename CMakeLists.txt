# This is the top level cmake
cmake_minimum_required(VERSION 3.28)


# ------------------
# PROJECT DEFINITION
# ------------------
# Deployment Target on MacOS
if(APPLE)
    # Force the deployment target to macOS 14
    set(CMAKE_OSX_DEPLOYMENT_TARGET "14.0" CACHE STRING "Minimum macOS version" FORCE)
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
option(APRIL_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(APRIL_BUILD_EXAMPLES   "Build examples" OFF)

# Backend Options (Bridges to C++ Preprocessor)
option(APRIL_ENABLE_XSIMD  "Use xsimd for vectorization" ON) # by default we use xsimd since std::simd is not mature yet on all platforms
option(APRIL_ENABLE_OPENMP "Use OpenMP executor backend" OFF)

# Code Health Options
option(APRIL_ENABLE_COVERAGE "Enable test coverage (GCC/Clang)" OFF)
option(APRIL_USE_ASAN        "Enable Address & Leak Sanitizer" OFF)
option(APRIL_USE_UBSAN       "Enable Undefined Behavior Sanitizer" OFF)
option(APRIL_USE_TSAN        "Enable Thread Sanitizer" OFF)

# Hardware & Optimization
option(APRIL_ENABLE_FAST_MATH "Enable unsafe floating-point optimizations (-ffast-math / /fp:fast)" ON)
set(APRIL_MARCH "native" CACHE STRING "Target architecture for GCC/Clang (e.g., native, haswell, x86-64-v3)")
set(APRIL_MSVC_ARCH "AVX2" CACHE STRING "Target instruction set for MSVC (e.g., AVX2, AVX512)")
set(APRIL_MSVC_FAVOR "NONE" CACHE STRING "Optimize for specific MSVC architecture (NONE, INTEL64, AMD64)")
set_property(CACHE APRIL_MSVC_FAVOR PROPERTY STRINGS "NONE" "INTEL64" "AMD64")

# Mutual Exclusion Checks
if(APRIL_USE_ASAN AND APRIL_USE_TSAN)
    message(FATAL_ERROR "[APRIL] AddressSanitizer (ASan) and ThreadSanitizer (TSan) are mutually exclusive. Pick one.")
endif()
if(APRIL_ENABLE_COVERAGE AND (APRIL_USE_ASAN OR APRIL_USE_TSAN OR APRIL_USE_UBSAN))
    message(WARNING "[APRIL] Coverage is enabled. Sanitizers will be ignored.")
endif()



# ------------------------
# APRIL LIBRARY DEFINITION
# ------------------------
# define header only library library with name April. This is what targets can link against
add_library(April INTERFACE)
target_compile_features(April INTERFACE cxx_std_23) # enforce c++23 on any target that uses April

# IGlob headers and attach them to the target for IDE
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

# Resolve SIMD Backend
if(APRIL_ENABLE_XSIMD)
    include(FetchContent) # get xSIMD
    FetchContent_Declare(
        xsimd
        GIT_REPOSITORY
        https://github.com/xtensor-stack/xsimd.git GIT_TAG 14.0.0
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
    target_compile_definitions(April INTERFACE APRIL_EXECUTOR_BACKEND_OMP)
endif()


# -----------
# BUILD FLAGS
# -----------
# This target holds all warnings, optimizations, and sanitizers.
# Internal targets (tests, examples) link to this. Users consuming the library do not.
add_library(April_dev_configs INTERFACE)

# Warnings
target_compile_options(April_dev_configs INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Werror -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W2 /permissive- /WX>
)

# Optimizations
target_compile_options(April_dev_configs INTERFACE
    # Release / RelWithDebInfo (GCC/Clang)
    $<$<AND:$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>,$<CXX_COMPILER_ID:GNU,Clang>>:
        -O3                      # Maximum optimization level
        -march=${APRIL_MARCH}    # Use all CPU instructions available on specified platform
        -fomit-frame-pointer     # Omit frame pointer for slightly faster calls
        -fno-plt                 # Avoid procedure linkage table, improves direct calls
        $<$<BOOL:${APRIL_ENABLE_FAST_MATH}>:-ffast-math> # Allow unsafe floating-point optimizations for speed
    >
    # Release / RelWithDebInfo (MSVC)
    $<$<AND:$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>,$<CXX_COMPILER_ID:MSVC>>:
        /O2                 # Maximize speed (replaces /Ox)
#        /GL                 # Whole program optimization (LTCG) # commented out because of long build times
        /Gw                 # Optimize global data
        /Ob3                # Aggressive inlining (better than /Ob2)
        /Zc:inline          # Remove unreferenced COMDATs (shrinks binary/link time)
        /cgthreads8         # Use 8 threads for compilation/codegen
        /GS-                # disable stack security checks: could lead to slight perf increase in hot loops
        /arch:${APRIL_MSVC_ARCH}  # Use all CPU instructions available on specified platform
        $<$<BOOL:${APRIL_ENABLE_FAST_MATH}>:/fp:fast> # Allow unsafe floating-point optimizations for speed
        $<$<NOT:$<STREQUAL:${APRIL_MSVC_FAVOR},NONE>>:/favor:${APRIL_MSVC_FAVOR}> # specify favoured architecture
        /bigobj
        /MP
    >

    # Debug
    $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:GNU,Clang>>:-O0 -g>
    $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:MSVC>>:/Od /Zi /bigobj>

    # Profile (GCC/Clang)
    $<$<AND:$<CONFIG:Profile>,$<CXX_COMPILER_ID:GNU,Clang>>:-O1 -g -pg>

    # Assembly (GCC/Clang)
    # Usage: Select "Assembly" build type, then run "make <target>.s"
    $<$<AND:$<CONFIG:Assembly>,$<CXX_COMPILER_ID:GNU,Clang>>:
        -O3
        -g                    # Generate debug info (for .loc annotations in asm)
        -fverbose-asm         # Add comments with variable names/source lines
        -masm=intel           # Enable Intel syntax
        -march=${APRIL_MARCH}
        -fomit-frame-pointer
        -finline-functions
        -fno-plt
        $<$<BOOL:${APRIL_ENABLE_FAST_MATH}>:-ffast-math>
    >
)

# Linker flags specifically for Profile mode
target_link_options(April_dev_configs INTERFACE
        $<$<AND:$<CONFIG:Profile>,$<CXX_COMPILER_ID:GNU,Clang>>:-pg>
)

# Sanitizers (Only if coverage is OFF)
if(NOT APRIL_ENABLE_COVERAGE)
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

if(APRIL_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
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
else()
    set(_active_arch "${APRIL_MARCH}")
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
message(STATUS "  Build Benchmarks:  ${APRIL_BUILD_BENCHMARKS}")
message(STATUS "  Build Examples:    ${APRIL_BUILD_EXAMPLES}")
message(STATUS "===========================================================")
