TEST(ComputeVEFTest, ClassIExample) {
    Vector3i result = PolyhedralLibrary::ComputeVEF(5, 1, 0);
    EXPECT_EQ(result[0], 12);  // V
    EXPECT_EQ(result[1], 30);  // E
    EXPECT_EQ(result[2], 20);  // F
}

TEST(ComputeVEFTest, ClassIIExample) {
    Vector3i result = PolyhedralLibrary::ComputeVEF(5, 2, 2);
    EXPECT_EQ(result[0], 3 * 4 * 20 + 2);  // V = 3*T*T*numF + 2, con T=4
    EXPECT_EQ(result[1], 9 * 4 * 20);      // E = 9*T*T*numF
    EXPECT_EQ(result[2], 6 * 4 * 20);      // F = 6*T*T*numF
}