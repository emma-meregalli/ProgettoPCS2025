#pragma once

#include <iostream>
#include <Eigen/Dense>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace PolyhedralLibrary
{

Eigen::Vector3i ComputeVEF(unsigned int q, int b, int c);

void CreateTxtFiles(const PolyhedralMesh& mesh);

bool GenerateDual(const PolyhedralMesh& mesh, PolyhedralMesh& dualMesh);

bool ExportTetrahedron(PolyhedralMesh& mesh, PolyhedralMesh& triMesh, const int& b, const int& c);

bool ExportOctahedron(PolyhedralMesh& mesh, PolyhedralMesh& triMesh, const int& b, const int& c);

bool ExportIcosahedron(PolyhedralMesh& mesh, PolyhedralMesh& triMesh, const int& b, const int& c);

void ExportParaView(PolyhedralMesh& mesh, bool cammino);

void OrderFaces(const vector<int>& unordered_faces, vector<int>& ordered_faces, const PolyhedralMesh& Mesh);

bool ShortestPath(PolyhedralMesh& mesh, unsigned int v1, unsigned int v2, unsigned int num_lati_iniziali, bool all_edges);

}
