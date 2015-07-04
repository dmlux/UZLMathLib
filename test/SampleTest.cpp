#include "gtest/gtest.h"
#include <uzlmath>

using namespace uzlmath;

TEST(SampleTest, AssertionTrue) {
    
    matrix<double> A(5, 5);
    rand(A, -1, 1);
    
    ASSERT_EQ(1, 1);
}