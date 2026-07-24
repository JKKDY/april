include(FetchContent)

# Fetch GoogleTest
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.16.0
    SYSTEM
)
set(gtest_build_shared_libs OFF CACHE BOOL "Build GTest as static" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(test_april
        test_april.cpp
        monitors/output_test.cpp
        monitors/xyz_test.cpp

        core/system_test.cpp
        core/regression_test.cpp

        containers/directsum_test.cpp
        containers/linkedcells_test.cpp
        containers/cell_ordering_test.cpp
        containers/scheduling_test.cpp

        exec/executors_test.cpp

        integrators/conservation_test.cpp
        integrators/stoermerverlet_test.cpp
        integrators/yoshida_test.cpp

        boundaries/apply_boundary_test.cpp
        boundaries/absorb_test.cpp
        boundaries/periodic_test.cpp
        boundaries/repulsive_test.cpp
        boundaries/reflective_test.cpp

        controllers/apply_controller_test.cpp
        controllers/thermostat_test.cpp

        fields/uniform_field_test.cpp
        fields/local_field_test.cpp
        fields/apply_field_test.cpp

        particle/particle_test.cpp
        particle/particle_fields_test.cpp
        particle/particle_generators_test.cpp
        particle/packed_fields_test.cpp

        math/range_test.cpp
        math/simd_test.cpp
        math/vec3_test.cpp
        math/simd_ref_test.cpp
        math/vec3_ref_test.cpp
        math/simd_vec3_ref_tests.cpp

        forces/interaction_test.cpp
        forces/forces_test.cpp
        forces/simd_forces_test.cpp

        utility/graph_test.cpp
)

target_include_directories(test_april PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/common)

# Link everything required
target_link_libraries(test_april PRIVATE
        April
        April_dev_configs
        GTest::gtest
        GTest::gtest_main
        GTest::gmock
        GTest::gmock_main
)

target_precompile_headers(test_april PRIVATE
        <gtest/gtest.h>
        <gmock/gmock.h>
        <vector>
        <memory>
)

# Discover GTest tests automatically
include(GoogleTest)
gtest_discover_tests(test_april)
set_target_properties(test_april PROPERTIES GTEST_DISCOVER_TESTS_MODE PRE_TEST)