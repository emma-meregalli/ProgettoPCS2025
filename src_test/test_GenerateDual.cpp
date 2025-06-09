#include <gtest/gtest.h>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <Eigen/Dense>

using namespace PolyhedralLibrary;
using namespace Eigen;

TEST(DualTest, DualMesh_Properties) {
    PolyhedralMesh mesh, triMesh, dual;
    ExportTetrahedron(mesh, triMesh, 1, 0);
    GenerateDual(triMesh, dual);
    EXPECT_EQ(dual.NumCell0Ds, triMesh.NumCell2Ds);
    EXPECT_EQ(dual.NumCell2Ds, triMesh.NumCell0Ds);
    EXPECT_GT(dual.NumCell1Ds, 0);
}