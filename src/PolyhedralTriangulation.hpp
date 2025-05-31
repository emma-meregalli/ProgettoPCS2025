#pragma once

#include <iostream>
#include "PolygonalMesh.hpp"

using namespace std;

namespace PolyhedralTriangulation
{

bool VertexIsDupe(const PolygonalMesh& mesh, const Vector3d& v);

bool EdgeIsDupe(const PolygonalMesh& mesh, const Vector2i& e);

void GenerateTriangulatedMesh(PolyhedralMesh& baseMesh, PolyhedralMesh& triMesh, unsigned int b, unsigned int c, const Vector3i& triDimensions);

}
