# top level 

cmake_minimum_required(VERSION 3.15)
project(April LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

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
            $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang>>:-O3 -DNDEBUG -march=native -ffast-math>
            # MSVC
            $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:
            /O2           # optimize for speed
            /GL           # whole program optimization (LTCG)
            /fp:fast      # fast-math (matches -ffast-math spirit)
            /arch:AVX2    # enable AVX2 on x64
            /favor:INTEL64
            /Oi           # intrinsic functions
            /Ot           # favor fast code
            /Ob3          # aggressive inlining
            /DNDEBUG
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
    if (NOT CMAKE_BUILD_TYPE MATCHES "Debug")
        message(STATUS "Skipping AddressSanitizer for target ${target_name} because build type is not Debug.")
        return()
    endif()

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        message(STATUS "Enabling AddressSanitizer + UBSan + Leak for target ${target_name} (GCC/Clang, Debug build)")
        target_compile_options(${target_name}
                PRIVATE
                -fsanitize=address,undefined,leak -fno-omit-frame-pointer
        )
        target_link_options(${target_name}
                PRIVATE
                -fsanitize=address,undefined,leak
        )
    elseif (MSVC)
        message(WARNING "Skipping AddressSanitizer on MSVC for target ${target_name} to avoid LNK2038 mismatches. Enable full sanitized build instead if desired.")
    endif()
endfunction()



add_compile_options(${GLOBAL_OPT_FLAGS})

# Include subdirectories
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(test)
add_subdirectory(benchmark)


