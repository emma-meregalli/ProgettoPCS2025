#include <gtest/gtest.h>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <Eigen/Dense>

using namespace PolyhedralLibrary;
using namespace Eigen;

TEST(ExportTest, ExportTetrahedronWorks) {
    PolyhedralMesh mesh, triMesh;
    EXPECT_TRUE(ExportTetrahedron(mesh, triMesh, 1, 0));
    EXPECT_GT(triMesh.NumCell0Ds, 0);
    EXPECT_GT(triMesh.NumCell1Ds, 0);
}

TEST(ExportTest, ExportOctahedronWorks) {
    PolyhedralMesh mesh, triMesh;
    EXPECT_TRUE(ExportOctahedron(mesh, triMesh, 1, 0));
    EXPECT_GT(triMesh.NumCell0Ds, 0);
    EXPECT_GT(triMesh.NumCell2Ds, 0);
}