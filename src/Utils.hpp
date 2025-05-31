#pragma once

#include <iostream>
#include "PolygonalMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{

std::Vector3i ComputeVEF(unsigned int q, int b, int c);

void CreateTxtFiles(const PolyhedralMesh& mesh);

bool GenerateDual(PolyhedralMesh& mesh, PolyhedralMesh& dualMesh);

bool ExportTetrahedron(PolyhedralMesh& mesh, PolyhedralMesh& triMesh, const int& b, const int& c);

bool ExportOctahedron(PolyhedralMesh& mesh, PolyhedralMesh& triMesh, const int& b, const int& c);

bool ExportIcosahedron(PolyhedralMesh& mesh, PolyhedralMesh& triMesh, const int& b, const int& c);

}
