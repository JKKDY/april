enable_testing()

# Use GTest from vcpkg
include(FetchContent)

# Download GoogleTest via FetchContent
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.16.0
)
FetchContent_MakeAvailable(googletest)



add_executable(test_april
    test_april.cpp
    environment_test.cpp
    map_test.cpp
)

target_link_libraries(test_april PRIVATE
    April
    GTest::gtest
    GTest::gtest_main
    GTest::gmock 
    GTest::gmock_main
)

include(GoogleTest)
gtest_discover_tests(test_april)
