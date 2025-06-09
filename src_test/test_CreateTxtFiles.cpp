#include <gtest/gtest.h>
#include "Utils.hpp"
#include <fstream>
#include <filesystem>

using namespace PolyhedralLibrary;
namespace fs = std::filesystem;

TEST(FilesTest, TxtFilesCreation) {
    PolyhedralMesh mesh, triMesh;
    ExportTetrahedron(mesh, triMesh, 1, 0);
    CreateTxtFiles(triMesh);

    EXPECT_TRUE(fs::exists("Cell0Ds.txt"));
    EXPECT_TRUE(fs::exists("Cell1Ds.txt"));
    EXPECT_TRUE(fs::exists("Cell2Ds.txt"));
    EXPECT_TRUE(fs::exists("Cell3Ds.txt"));
}