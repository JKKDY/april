enable_testing()

# Use GTest from vcpkg
find_package(GTest CONFIG REQUIRED)

add_executable(test_april
    test_april.cpp
)

target_link_libraries(test_april PRIVATE
    APRIL
    GTest::gtest
    GTest::gtest_main
    GTest::gmock 
    GTest::gmock_main
)

include(GoogleTest)
gtest_discover_tests(test_april)
