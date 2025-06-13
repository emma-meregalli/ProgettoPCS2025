#include <gtest/gtest.h>
#include "PolyhedralTriangulation.hpp"
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include <Eigen/Dense>

using namespace PolyhedralTriangulation;
using namespace PolyhedralLibrary;
using namespace Eigen;

TEST(TriangulationTest, GenerateTriangulatedMesh2_ClassII) {
    PolyhedralMesh baseMesh, triMesh;

    // Definizione di un triangolo di base
    baseMesh.Cell0DsCoordinates = MatrixXd(3, 3);
    baseMesh.Cell0DsCoordinates << 0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0;
    baseMesh.Cell0DsId = {0, 1, 2};

    baseMesh.Cell2DsId = {0};
    baseMesh.Cell2DsVertices = {{0, 1, 2}};

    // Simula la dimensione della mesh triangolata
    Vector3i VEF = PolyhedralLibrary::ComputeVEF(3, 1, 1);  // Classe II: b == c

    bool success = GenerateTriangulatedMesh2(baseMesh, triMesh, 1, 1, VEF);

    EXPECT_TRUE(success);
    EXPECT_GT(triMesh.Cell0DsId.size(), 3);  // Vertici triangolati
    EXPECT_GT(triMesh.Cell1DsId.size(), 0);  // Lati triangolati
    EXPECT_GT(triMesh.Cell2DsId.size(), 0);  // Facce triangolate
}