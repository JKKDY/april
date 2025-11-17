#enable_testing()

include(FetchContent)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.16.0
)

set(gtest_build_shared_libs OFF CACHE BOOL "Build GTest as static")
FetchContent_MakeAvailable(googletest)

add_executable(test_april
        test_april.cpp
        system_test.cpp
        forces/interaction_test.cpp
        containers/directsum_test.cpp
        output_test.cpp
        integrators/stoermerverlet_test.cpp
        containers/linkedcells_test.cpp
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
        forces/forces_test.cpp
)

target_include_directories(test_april PRIVATE
        ${PROJECT_SOURCE_DIR}/test/common
)

target_link_libraries(test_april PRIVATE
        April
        GTest::gtest
        GTest::gtest_main
        GTest::gmock
        GTest::gmock_main
)

set_target_build_flags(test_april)
enable_sanitizer_if_dbg(test_april)
enable_coverage(test_april)


# Discover GTest tests
include(GoogleTest)
gtest_discover_tests(test_april)
set_target_properties(test_april PROPERTIES GTEST_DISCOVER_TESTS_MODE PRE_TEST)
