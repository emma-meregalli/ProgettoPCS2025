#pragma once

#include <iostream>
#include "PolygonalMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{

std::vector<int> ComputeVEF(unsigned int q, int b, int c);

void CreateTxtFiles(const PolyhedralMesh& mesh);

bool ExportDual(PolyhedralMesh& mesh, PolyhedralMesh& dualMesh);

bool ExportTetrahedron(PolyhedralMesh& mesh);

bool ExportCube(PolyhedralMesh& mesh);

bool ExportOctahedron(PolyhedralMesh& mesh);

bool ExportDodecahedron(PolyhedralMesh& mesh);

bool ExportIcosahedron(PolyhedralMesh& mesh);

}
