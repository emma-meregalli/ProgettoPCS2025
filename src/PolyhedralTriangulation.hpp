#pragma once

#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace PolyhedralLibrary;

namespace PolyhedralTriangulation
{

bool VertexIsDupe(const PolyhedralMesh& mesh, const Vector3d& v, unsigned int& original_id);

bool EdgeIsDupe(const PolyhedralMesh& mesh, const Vector2i& e);

bool GenerateTriangulatedMesh1(PolyhedralMesh& baseMesh, PolyhedralMesh& triMesh, const unsigned int& b, const unsigned int& c, const Vector3i& triDimensions);

bool GenerateTriangulatedMesh2(PolyhedralMesh& baseMesh, PolyhedralMesh& triMesh, const unsigned int& b, const unsigned int& c, const Vector3i& triDimensions);

}
