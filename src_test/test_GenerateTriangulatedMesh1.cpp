#include <gtest/gtest.h>
#include "PolyhedralTriangulation.hpp"
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include <Eigen/Dense>

using namespace PolyhedralTriangulation;
using namespace PolyhedralLibrary;
using namespace Eigen;

TEST(TriangulationTest, GenerateTriangulatedMesh1_BasicTriangle) {
    PolyhedralMesh baseMesh, triMesh;

    // Vertici di un triangolo semplice
    baseMesh.Cell0DsCoordinates = MatrixXd(3, 3);
    baseMesh.Cell0DsCoordinates << 0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0;
    baseMesh.Cell0DsId = {0, 1, 2};
    
    // Una faccia triangolare
    baseMesh.Cell2DsId = {0};
    baseMesh.Cell2DsVertices = {{0, 1, 2}};

    // Simula il VEF da ComputeVEF
    Vector3i VEF = PolyhedralLibrary::ComputeVEF(3, 1, 0);

    bool success = GenerateTriangulatedMesh1(baseMesh, triMesh, 1, 0, VEF);

    EXPECT_TRUE(success);
    EXPECT_GT(triMesh.Cell0DsId.size(), 3);  // Devono essere generati nuovi vertici
    EXPECT_GT(triMesh.Cell1DsId.size(), 0);  // Devono essere generati lati
    EXPECT_GT(triMesh.Cell2DsId.size(), 0);  // Devono essere generate facce
}
