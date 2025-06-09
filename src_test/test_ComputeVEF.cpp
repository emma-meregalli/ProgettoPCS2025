#include <gtest/gtest.h>
#include "Utils.hpp"
#include <Eigen/Dense>

using namespace Eigen;
using namespace PolyhedralLibrary;

TEST(UtilsTest, ComputeVEF_Tetrahedron) {
    Vector3i result = ComputeVEF(3, 1, 0);
    EXPECT_EQ(result[0], 4);
    EXPECT_EQ(result[1], 6);
    EXPECT_EQ(result[2], 4);
}

TEST(UtilsTest, ComputeVEF_ZeroSubdivision) {
    Vector3i result = ComputeVEF(3, 0, 0);
    EXPECT_EQ(result[0], 2);
    EXPECT_EQ(result[1], 0);
    EXPECT_EQ(result[2], 0);
}

TEST(UtilsTest, ComputeVEF_Icosahedron) {
    Vector3i result = ComputeVEF(5, 2, 2);
    EXPECT_EQ(result[0], 122); // 10*T + 2, T=12
    EXPECT_EQ(result[1], 360); // 30*T
    EXPECT_EQ(result[2], 240); // 20*T
}