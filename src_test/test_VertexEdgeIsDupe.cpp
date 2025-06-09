#include <gtest/gtest.h>
#include "PolyhedralTriangulation.hpp"
#include "PolyhedralMesh.hpp"
#include <Eigen/Dense>

using namespace PolyhedralTriangulation;
using namespace PolyhedralLibrary;
using namespace Eigen;

TEST(TriangulationTest, VertexIsDupe_True) {
    PolyhedralMesh mesh;
    mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 1);
    mesh.Cell0DsCoordinates.col(0) << 1.0, 0.0, 0.0;
    mesh.Cell0DsId = {0};
    Vector3d v(1.0 + 1e-13, 0.0, 0.0);
    EXPECT_TRUE(VertexIsDupe(mesh, v));
}

TEST(TriangulationTest, VertexIsDupe_False) {
    PolyhedralMesh mesh;
    mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 1);
    mesh.Cell0DsCoordinates.col(0) << 1.0, 0.0, 0.0;
    mesh.Cell0DsId = {0};
    Vector3d v(0.0, 1.0, 0.0);
    EXPECT_FALSE(VertexIsDupe(mesh, v));
}

TEST(TriangulationTest, EdgeIsDupe_True) {
    PolyhedralMesh mesh;
    mesh.Cell1DsId = {0};
    mesh.Cell1DsExtrema = MatrixXi(2, 1);
    mesh.Cell1DsExtrema << 1, 2;
    Vector2i edge(2, 1);
    EXPECT_TRUE(EdgeIsDupe(mesh, edge));
}

TEST(TriangulationTest, EdgeIsDupe_False) {
    PolyhedralMesh mesh;
    mesh.Cell1DsId = {0};
    mesh.Cell1DsExtrema = MatrixXi(2, 1);
    mesh.Cell1DsExtrema << 0, 1;
    Vector2i edge(2, 3);
    EXPECT_FALSE(EdgeIsDupe(mesh, edge));
}