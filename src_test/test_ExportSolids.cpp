#include <gtest/gtest.h>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <Eigen/Dense>

using namespace PolyhedralLibrary;
using namespace Eigen;

TEST(ExportTetrahedronTest, BasicExportWorks) {
    PolyhedralMesh mesh, triMesh;
    bool success = PolyhedralLibrary::ExportTetrahedron(mesh, triMesh, 1, 0);
    EXPECT_TRUE(success);
    EXPECT_EQ(mesh.Cell0DsId.size(), 4);
    EXPECT_EQ(mesh.Cell1DsId.size(), 6);
    EXPECT_EQ(mesh.Cell2DsId.size(), 4);
}

TEST(ExportSolidsTest, ExportOctahedronProducesCorrectStructure) {
    PolyhedralMesh mesh, triMesh;
    bool success = PolyhedralLibrary::ExportOctahedron(mesh, triMesh, 1, 0);
    EXPECT_TRUE(success);
    EXPECT_EQ(mesh.Cell0DsId.size(), 6);
    EXPECT_EQ(mesh.Cell1DsId.size(), 12);
    EXPECT_EQ(mesh.Cell2DsId.size(), 8);
    EXPECT_EQ(mesh.Cell3DsId.size(), 1);
}

TEST(ExportSolidsTest, ExportIcosahedronProducesCorrectStructure) {
    PolyhedralMesh mesh, triMesh;
    bool success = PolyhedralLibrary::ExportIcosahedron(mesh, triMesh, 1, 0);
    EXPECT_TRUE(success);
    EXPECT_EQ(mesh.Cell0DsId.size(), 12);
    EXPECT_EQ(mesh.Cell1DsId.size(), 30);
    EXPECT_EQ(mesh.Cell2DsId.size(), 20);
    EXPECT_EQ(mesh.Cell3DsId.size(), 1);
}
