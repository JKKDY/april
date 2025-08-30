enable_testing()

# Use GTest from vcpkg
include(FetchContent)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.16.0
)
FetchContent_MakeAvailable(googletest)

add_executable(test_april
        test_april.cpp
        system_test.cpp
        set_test.cpp
        interaction_test.cpp
        directsum_test.cpp
        output_test.cpp
        stoermerverlet_test.cpp
        linkedcells_test.cpp
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

# Discover GTest tests
include(GoogleTest)
gtest_discover_tests(test_april)
