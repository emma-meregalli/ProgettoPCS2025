#pragma once

#include <iostream>
#include "PolygonalMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{

bool ExportTetrahedron(PolyhedralMesh& mesh);

bool ExportCube(PolyhedralMesh& mesh);

bool ExportOctahedron(PolyhedralMesh& mesh);

bool ExportDodecahedron(PolyhedralMesh& mesh);

bool ExportIcosahedron(PolyhedralMesh& mesh);

}
