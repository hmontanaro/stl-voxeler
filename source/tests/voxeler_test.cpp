#include <gtest/gtest.h>
#include "../voxeler.h" 

// Add Vector3 tests
TEST(Vector3Test, DefaultConstructor) {
    Vector3 v;
    EXPECT_DOUBLE_EQ(0.0, v.x);
    EXPECT_DOUBLE_EQ(0.0, v.y);
    EXPECT_DOUBLE_EQ(0.0, v.z);
}

TEST(Vector3Test, Addition) {
    Vector3 v1(1.0, 2.0, 3.0);
    Vector3 v2(4.0, 5.0, 6.0);
    Vector3 result = v1 + v2;

    EXPECT_DOUBLE_EQ(5.0, result.x);
    EXPECT_DOUBLE_EQ(7.0, result.y);
    EXPECT_DOUBLE_EQ(9.0, result.z);
}

// Define a simple test fixture for triangle-box tests
class TriBoxTest : public ::testing::Test {
protected:
    // Create a unit box centered at origin
    Vector3 boxCenter{ 0, 0, 0 };
    Vector3 boxHalfSize{ 0.5, 0.5, 0.5 };

    // Set up any other common test data
};

TEST_F(TriBoxTest, TriangleInsideBox) {
    Triangle tri{
        Vector3(-0.1, -0.1, -0.1),
        Vector3(0.1,  0.1, -0.1),
        Vector3(0.0,  0.0,  0.1)
    };

    EXPECT_TRUE(triBoxOverlap(boxCenter, boxHalfSize, tri));
}

TEST_F(TriBoxTest, TriangleOutsideBox) {
    Triangle tri{
        Vector3(1.0, 1.0, 1.0),
        Vector3(1.1, 1.1, 1.1),
        Vector3(1.2, 1.0, 1.0)
    };

    EXPECT_FALSE(triBoxOverlap(boxCenter, boxHalfSize, tri));
}

TEST_F(TriBoxTest, TriangleIntersectingBox) {
    Triangle tri{
        Vector3(-1.0, -1.0, 0.0),
        Vector3(1.0, -1.0, 0.0),
        Vector3(0.0,  1.0, 0.0)
    };

    EXPECT_TRUE(triBoxOverlap(boxCenter, boxHalfSize, tri));
}
