#include <gtest/gtest.h>
#include <april/april.h>

TEST(APRILTest, HelloRuns) {
    EXPECT_NO_THROW(april::hello());
}