#pragma once

#include <iostream>
#include "PolygonalMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{

bool ImportMesh(PolyhedralMesh& mesh);

bool ImportCell0Ds(PolyhedralMesh& mesh);

bool ImportCell1Ds(PolyhedralMesh& mesh);

bool ImportCell2Ds(PolyhedralMesh& mesh);

bool CheckMarkers(const PolyhedralMesh& mesh);

bool CheckEdges(const PolyhedralMesh& mesh);

bool CheckAreas(const PolyhedralMesh& mesh);

}
