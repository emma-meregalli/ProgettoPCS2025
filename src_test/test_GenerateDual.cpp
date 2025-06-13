#include <gtest/gtest.h>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <Eigen/Dense>

using namespace PolyhedralLibrary;
using namespace Eigen;

TEST(GenerateDualTest, GeneratesDualCorrectly) {
    PolyhedralMesh mesh, dual;
    PolyhedralLibrary::ExportTetrahedron(mesh, mesh, 1, 0);
    bool generated = PolyhedralLibrary::GenerateDual(mesh, dual);
    EXPECT_TRUE(generated);
    EXPECT_EQ(dual.NumCell0Ds, mesh.Cell2DsId.size());
    EXPECT_EQ(dual.NumCell2Ds, mesh.Cell0DsId.size());
}