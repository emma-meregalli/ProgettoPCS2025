#pragma once

#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace PolyhedralTriangulation
{

bool VertexIsDupe(const PolygonalMesh& mesh, const Vector3d& v);

bool EdgeIsDupe(const PolygonalMesh& mesh, const Vector2i& e);

bool GenerateTriangulatedMesh(PolyhedralMesh& baseMesh, PolyhedralMesh& triMesh, const unsigned int& b, const unsigned int& c, const Vector3i& triDimensions);

}
