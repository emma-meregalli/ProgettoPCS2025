TEST(CreateTxtFilesTest, GeneratesAllTxtFiles) {
    PolyhedralMesh mesh;
    
    mesh.Cell0DsId = {0};
    mesh.Cell0DsCoordinates = MatrixXd::Zero(3,1);
    mesh.Cell1DsId = {0};
    mesh.Cell1DsExtrema = MatrixXi::Zero(2,1);
    mesh.Cell2DsId = {0};
    mesh.Cell2DsVertices = {{0}};
    mesh.Cell2DsEdges = {{0}};
    mesh.Cell3DsId = {0};
    mesh.Cell3DsVertices = {0};
    mesh.Cell3DsEdges = {0};
    mesh.Cell3DsFaces = {0};

    PolyhedralLibrary::CreateTxtFiles(mesh);

    std::ifstream file("Cell0Ds.txt");
    ASSERT_TRUE(file.is_open());
    std::string line;
    std::getline(file, line);
    EXPECT_EQ(line, "ID;x;y;z");
}
