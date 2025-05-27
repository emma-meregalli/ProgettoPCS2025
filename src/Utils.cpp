#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
vector<int> ComputeVEF(unsigned int q, int b, int c)
{
    vector<int> VEF(3, 0);  // inizializzo un vettore nullo di lunghezza 3

    unsigned int T = b * b + b * c + c * c;
    unsigned int V, E, F;

    if (q == 3) {
        V = 2 * T + 2;
        E = 6 * T;
        F = 4 * T;
    } else if (q == 4) {
        V = 4 * T + 2;
        E = 12 * T;
        F = 8 * T;
    } else {
        V = 10 * T + 2;
        E = 30 * T;
        F = 20 * T;
    }

    VEF[0] = V;  // V = vertici
    VEF[1] = E;  // E = spigoli
    VEF[2] = F;  // F = facce

    return VEF;  // Mi viene restituito il vettore con i valori di V, E, F
}

void CreateTxtFiles(const PolyhedralMesh& mesh) {
    // Creazione Cell0Ds.txt
    ofstream Cell0Ds("Cell0Ds.txt");
    out0 << "ID;x;y;z\n";
    for (size_t i = 0; i < mesh.Cell0DsId.size(); i++) {
        Cell0Ds << mesh.Cell0DsId[i] << ";" << mesh.Cell0DsCoordinates(0, i) << ";" << mesh.Cell0DsCoordinates(1, i) << ";" << mesh.Cell0DsCoordinates(2, i) << "\n";
    }
    Cell0Ds.close();

    // Creazione Cell1Ds.txt 
    ofstream Cell1Ds("Cell1Ds.txt");
    Cell1Ds << "ID;Origin;End\n";
    for (size_t i = 0; i < mesh.Cell1DsId.size(); i++) {
        Cell1Ds << mesh.Cell1DsId[i] << ";"
             << mesh.Cell1DsExtrema(0, i) << ";"
             << mesh.Cell1DsExtrema(1, i) << "\n";
    }
    Cell1Ds.close();

    // Creazione Cell2Ds.txt 
    ofstream Cell2Ds("Cell2Ds.txt");
    Cell2Ds << "ID;NumVertices;NumEdges;Vertices;Edges\n";
    for (size_t i = 0; i < mesh.Cell2DsId.size(); i++) {
        Cell2Ds << mesh.Cell2DsId[i] << ";"
             << mesh.Cell2DsVertices[i].size() << ";"
             << mesh.Cell2DsEdges[i].size();

        // Vertici
        for (unsigned int v : mesh.Cell2DsVertices[i])
            Cell2Ds << ";" << v;

        // Lati
        for (unsigned int e : mesh.Cell2DsEdges[i])
            Cell2Ds << ";" << e;

        Cell2Ds << "\n";
    }
    Cell2Ds.close();

    // Creazione Cell3Ds.txt
    ofstream Cell3Ds("Cell3Ds.txt");
    Cell3Ds << "ID;Vertices;Edges;Faces\n";
    for (size_t i = 0; i < mesh.Cell3DsId.size(); i++) {
        Cell3Ds << mesh.Cell3DsId[i];

        // Vertici
        for (unsigned int v : mesh.Cell3DsVertices[i])
            Cell3Ds << ";" << v;

        // Lati
        for (unsigned int e : mesh.Cell3DsEdges[i])
            Cell3Ds << ";" << e;

        // Facce
        for (unsigned int f : mesh.Cell3DsFaces[i])
            Cell3Ds << ";" << f;

        Cell3Ds << "\n";
    }
    Cell3Ds.close();
} 

bool ExportDual(PolyhedralMesh& mesh, PolyhedralMesh& dualMesh){
    MatrixXd barycenters = {};
    barycenters.reserve(mesh.NumCell2Ds,3);
    Vector3d tmp = zeros();
    for (unsigned int i=0; i<mesh.NumCell2Ds; i++){
        for (unsigned int v : mesh.Cell2DsVertices[i]){
            tmp += mesh.Cell0DsCoordinates.row(v);
        }
        tmp /= mesh.Cell2DsVertices[i].size(); 
        barycenters.push_back(tmp);
    }

    // Riempie Cell0Ds del duale con i baricentri
    int numDualVertices = barycenters.size();
    dualMesh.Cell0DsCoordinates.resize(numDualVertices,3);
    for (int i = 0; i < numDualVertices; ++i) {
        dualMesh.Cell0DsCoordinates.row(i) = barycenters[i];
        dualMesh.Cell0DsId.push_back(i);
    }

    // Mappa: vertice originale -> facce adiacenti che lo contengono
    map<int, vector<int>> vertexToFaces;
    for (size_t f = 0; f < mesh.Cell2DsVertices.size(); ++f) {
        for (unsigned int v : mesh.Cell2DsVertices[f]) {
            vertexToFaces[v].push_back(f);
        }
    }

    // Step 2: crea facce duali (una per ogni vertice originale)
    int faceId = 0;
    int edgeId = 0;
    map<pair<int, int>, int> edgeMap; // per evitare spigoli duplicati
    for (const auto& [vertex, faces] : vertexToFaces) {
        // Ordina ciclicamente le facce attorno al vertice originale
        vector<int> orderedFaces;
        vector<int> rest = faces;

        orderedFaces.push_back(rest[0]); // inizia da una qualsiasi
        rest.erase(rest.begin());

        while (!rest.empty()) {
            int current = orderedFaces.back();
            int next = -1;
            for (auto it = rest.begin(); it != rest.end(); ++it) {
                const auto& fv1 = mesh.Cell2DsVertices[current];
                const auto& fv2 = mesh.Cell2DsVertices[*it];
                vector<int> common;
                for (int a : fv1)
                    for (int b : fv2)
                        if (a == b) common.push_back(a);
                if (common.size() >= 2) {
                    next = *it;
                    rest.erase(it);
                    break;
                }
            }
            if (next == -1) break;
            orderedFaces.push_back(next);
        }

        // Crea la nuova faccia nel duale
        dualMesh.Cell2DsVertices.push_back(orderedFaces);
        dualMesh.Cell2DsId.push_back(faceId);

        // Costruzione degli spigoli per la faccia
        vector<unsigned int> faceEdges;
        for (size_t i = 0; i < orderedFaces.size(); ++i) {
            int a = orderedFaces[i];
            int b = orderedFaces[(i + 1) % orderedFaces.size()];
            pair<int, int> key = (a < b) ? make_pair(a, b) : make_pair(b, a);
            if (edgeMap.count(key) == 0) {
                dualMesh.Cell1DsExtrema.conservativeResize(2, edgeId + 1);
                dualMesh.Cell1DsExtrema(0, edgeId) = key.first;
                dualMesh.Cell1DsExtrema(1, edgeId) = key.second;
                dualMesh.Cell1DsId.push_back(edgeId);
                edgeMap[key] = edgeId++;
            }
            faceEdges.push_back(edgeMap[key]);
        }
        dualMesh.Cell2DsEdges.push_back(faceEdges);
        ++faceId;
    }

    // Definisce il poliedro 3D che contiene tutto
    dualMesh.Cell3DsId = {0};
    dualMesh.Cell3DsVertices = dualMesh.Cell0DsId;
    dualMesh.Cell3DsEdges = dualMesh.Cell1DsId;
    dualMesh.Cell3DsFaces = dualMesh.Cell2DsId;

    // Aggiorna i conteggi delle celle
    dualMesh.NumCell0Ds = dualMesh.Cell0DsId.size();
    dualMesh.NumCell1Ds = dualMesh.Cell1DsId.size();
    dualMesh.NumCell2Ds = dualMesh.Cell2DsId.size();
    dualMesh.NumCell3Ds = 1;
}

bool ExportTetrahedron(PolyhedralMesh& mesh) {
	
	// Vertici
    double r = sqrt(3.0) / 3.0;
    
    mesh.Cell0DsCoordinates.reserve(4,3);
    mesh.Cell0DsId.reserve(4);

    mesh.Cell0DsCoordinates(0, 0) = r;   mesh.Cell0DsCoordinates(1, 0) = r;   mesh.Cell0DsCoordinates(2, 0) = r;  
    mesh.Cell0DsCoordinates(0, 1) = -r;  mesh.Cell0DsCoordinates(1, 1) = -r;  mesh.Cell0DsCoordinates(2, 1) = r; 
    mesh.Cell0DsCoordinates(0, 2) = -r;  mesh.Cell0DsCoordinates(1, 2) = r;   mesh.Cell0DsCoordinates(2, 2) = -r;
    mesh.Cell0DsCoordinates(0, 3) = r;   mesh.Cell0DsCoordinates(1, 3) = -r;  mesh.Cell0DsCoordinates(2, 3) = -r;
    
    mesh.Cell0DsId = {0,1,2,3};
    
    // Lati
    mesh.Cell1DsId.reserve(6);
    Cell1DsId = {0,1,2,3,4,5};

    mesh.Cell1DsExtrema.resize(6, 2);
    mesh.Cell1DsExtrema <<
        0, 1,
        1, 2,
        2, 0,
        0, 3,
        3, 1,
        2, 3;

    // Facce
    mesh.Cell2DsId.reserve(4);
    mesh.Cell2DsId = {0, 1, 2, 3};
    
    mesh.Cell2DsVertices.resize(4);
	for (auto& v : mesh.Cell2DsVertices)
		v.reserve(3);
	
	mesh.Cell2DsEdges.resize(4);
	for (auto& e : mesh.Cell2DsEdges)
		e.reserve(3);

    mesh.Cell2DsVertices = {
        {0, 1, 2}, 
        {0, 3, 1}, 
        {1, 3, 2},
        {2, 3, 0} 
    };

    mesh.Cell2DsEdges = {
        {0, 1, 2}, 
        {3, 4, 0},
        {4, 5, 1},
        {5, 3, 2} 
    };

    // Poliedro
    mesh.Cell3DsId.reserve(1);
    mesh.Cell3DsVertices.reserve(4);
    mesh.Cell3DsEdges.reserve(6);
    mesh.Cell3DsFaces.reserve(4);

    mesh.Cell3DsId = {0};
    mesh.Cell3DsVertices = {0, 1, 2, 3};
    mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5};
    mesh.Cell3DsFaces = {0, 1, 2, 3};

    CreateTxtFiles(const PolyhedralMesh& mesh);
    ExportDual(const PolyhedralMesh& mesh, PolyhedralMesh& dualMesh);
    CreateTxtFiles(const PolyhedralMesh& dualMesh);

    return true;
}

bool ExportCube(PolyhedralMesh& mesh) {
	
	// Vertici
	double r = sqrt(3.0) / 3.0;

	mesh.Cell0DsCoordinates.reserve(8, 3); 
	mesh.Cell0DsId.reserve(8);

	mesh.Cell0DsCoordinates(0, 0) = r;   mesh.Cell0DsCoordinates(1, 0) = r;   mesh.Cell0DsCoordinates(2, 0) = r;  
	mesh.Cell0DsCoordinates(0, 1) = r;   mesh.Cell0DsCoordinates(1, 1) = r;   mesh.Cell0DsCoordinates(2, 1) = -r;  
	mesh.Cell0DsCoordinates(0, 2) = r;   mesh.Cell0DsCoordinates(1, 2) = -r;  mesh.Cell0DsCoordinates(2, 2) = r;  
	mesh.Cell0DsCoordinates(0, 3) = r;   mesh.Cell0DsCoordinates(1, 3) = -r;  mesh.Cell0DsCoordinates(2, 3) = -r;  
	mesh.Cell0DsCoordinates(0, 4) = -r;  mesh.Cell0DsCoordinates(1, 4) = r;   mesh.Cell0DsCoordinates(2, 4) = r;  
	mesh.Cell0DsCoordinates(0, 5) = -r;  mesh.Cell0DsCoordinates(1, 5) = r;   mesh.Cell0DsCoordinates(2, 5) = -r;  
	mesh.Cell0DsCoordinates(0, 6) = -r;  mesh.Cell0DsCoordinates(1, 6) = -r;  mesh.Cell0DsCoordinates(2, 6) = r;  
	mesh.Cell0DsCoordinates(0, 7) = -r;  mesh.Cell0DsCoordinates(1, 7) = -r;  mesh.Cell0DsCoordinates(2, 7) = -r;
	
	mesh.Cell0DsId = {0,1,2,3,4,5,6,7};
	
	// Lati
    mesh.Cell1DsId.reserve(12);
    Cell1DsId = {0,1,2,3,4,5,6,7,8,9,10,11};

    mesh.Cell1DsExtrema = MatrixXi::Zero(12, 2);
	mesh.Cell1DsExtrema <<
		0, 1,  // Lato 0
		1, 2,  // Lato 1
		2, 3,  // Lato 2
		3, 0,  // Lato 3
		4, 5,  // Lato 4
		5, 6,  // Lato 5
		6, 7,  // Lato 6
		7, 4,  // Lato 7
		0, 4,  // Lato 8
		1, 5,  // Lato 9
		2, 6,  // Lato 10
		3, 7;  // Lato 11
		
	// Facce
	mesh.Cell2DsId.reserve(6); 
	mesh.Cell2DsId = {0, 1, 2, 3, 4, 5};
	
	mesh.Cell2DsVertices.resize(6);
	for (auto& v : mesh.Cell2DsVertices)
		v.reserve(4);
		
	mesh.Cell2DsEdges.resize(6);
	for (auto& e : mesh.Cell2DsEdges)
		e.reserve(4);
	
	mesh.Cell2DsVertices = {
		{0, 1, 2, 3},  
		{4, 5, 6, 7},  
		{0, 1, 5, 4},  
		{1, 2, 6, 5},  
		{2, 3, 7, 6},
		{3, 0, 4, 7} 
    };
    
	mesh.Cell2DsEdges = {
		{0, 1, 2, 3},
		{4, 5, 6, 7},
		{0, 9, 4, 8},
		{1, 10, 5, 9},
		{2, 11, 6, 10},
		{3, 8, 7, 11}
	};
	
	// Poliedro 
	mesh.Cell3DsId.reserve(1);
	mesh.Cell3DsVertices.reserve(8); 
	mesh.Cell3DsEdges.reserve(12);
	mesh.Cell3DsFaces.reserve(6); 
	
	mesh.Cell3DsId = {0};
	mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5, 6, 7};
	mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
	mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5};
	
    ExportDual(const PolyhedralMesh& mesh, PolyhedralMesh& dualMesh);
    CreateTxtFiles(const PolyhedralMesh& dualMesh);

	return true;
}
	
bool ExportOctahedron(PolyhedralMesh& mesh) {
	
    // Vertici
    double r = 1.0; 

    mesh.Cell0DsCoordinates = MatrixXd::Zero(6, 3);
    mesh.Cell0DsId.reserve(6);

    mesh.Cell0DsCoordinates(0, 0) = r;   mesh.Cell0DsCoordinates(1, 0) = 0.0;  mesh.Cell0DsCoordinates(2, 0) = 0.0;  
    mesh.Cell0DsCoordinates(0, 1) = -r;  mesh.Cell0DsCoordinates(1, 1) = 0.0;  mesh.Cell0DsCoordinates(2, 1) = 0.0;  
    mesh.Cell0DsCoordinates(0, 2) = 0.0;  mesh.Cell0DsCoordinates(1, 2) = r;   mesh.Cell0DsCoordinates(2, 2) = 0.0;  
    mesh.Cell0DsCoordinates(0, 3) = 0.0;  mesh.Cell0DsCoordinates(1, 3) = -r;  mesh.Cell0DsCoordinates(2, 3) = 0.0;  
    mesh.Cell0DsCoordinates(0, 4) = 0.0;  mesh.Cell0DsCoordinates(1, 4) = 0.0;  mesh.Cell0DsCoordinates(2, 4) = r;   
    mesh.Cell0DsCoordinates(0, 5) = 0.0;  mesh.Cell0DsCoordinates(1, 5) = 0.0;  mesh.Cell0DsCoordinates(2, 5) = -r;  

    mesh.Cell0DsId = {0, 1, 2, 3, 4, 5};

    // Lati
    mesh.Cell1DsId.reserve(12);
    Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    mesh.Cell1DsExtrema = MatrixXi::Zero(12, 2);
    mesh.Cell1DsExtrema <<
        0, 2,  // Lato 0
		2, 1,  // Lato 1
		1, 3,  // Lato 2
		3, 0,  // Lato 3
		0, 4,  // Lato 4
		2, 4,  // Lato 5
		4, 1,  // Lato 6
		4, 3,  // Lato 7
		1, 5,  // Lato 8
		3, 5,  // Lato 9
		0, 5,  // Lato 10
		2, 5;  // Lato 11
	
    // Facce
    mesh.Cell2DsId.reserve(8);
    mesh.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7};

    mesh.Cell2DsVertices.resize(8);
	for (auto& v : mesh.Cell2DsVertices)
		v.reserve(3);
		
	mesh.Cell2DsEdges.resize(8);
	for (auto& e : mesh.Cell2DsEdges)
		e.reserve(3);
	

    mesh.Cell2DsVertices = {
        {0, 2, 4},  // Faccia 0
		{0, 4, 3},  // Faccia 1
		{3, 4, 1},  // Faccia 2
		{1, 2, 4},  // Faccia 3
		{2, 5, 0},  // Faccia 4
		{2, 1, 5},  // Faccia 5
		{1, 5, 3},  // Faccia 6
		{0, 3, 5}   // Faccia 7
    };

    mesh.Cell2DsEdges = {
        {0, 5, 4},  // Faccia 0
		{4, 7, 3},  // Faccia 1
		{7, 6, 2}, // Faccia 2
		{1, 5, 6},  // Faccia 3
		{11, 10, 0}, // Faccia 4
		{1, 8 , 11}, // Faccia 5
		{8, 9, 2},  // Faccia 6
		{3, 9, 10}  // Faccia 7
    };

    // Poliedro
    mesh.Cell3DsId.reserve(1);
    mesh.Cell3DsVertices.reserve(6);
    mesh.Cell3DsEdges.reserve(12);
    mesh.Cell3DsFaces.reserve(8);

    mesh.Cell3DsId = {0};
    mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5};
    mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7};
    
    CreateTxtFiles(const PolyhedralMesh& mesh);

    return true;
}

bool ExportDodecahedron(PolyhedralMesh& mesh) {
	
    // Vertici
    double r = 1.0;
    double phi = (1.0 + sqrt(5.0)) / 2.0;

    mesh.Cell0DsCoordinates = MatrixXd::Zero(20, 3);
    mesh.Cell0DsId.reserve(20);
    double norm = sqrt(3.0);
    
    mesh.Cell0DsCoordinates(0, 0) = r / norm;  mesh.Cell0DsCoordinates(1, 0) = r / norm;  mesh.Cell0DsCoordinates(2, 0) = r / norm;
	mesh.Cell0DsCoordinates(0, 1) = r / norm;  mesh.Cell0DsCoordinates(1, 1) = r / norm;  mesh.Cell0DsCoordinates(2, 1) = -r / norm;
	mesh.Cell0DsCoordinates(0, 2) = r / norm;  mesh.Cell0DsCoordinates(1, 2) = -r / norm;  mesh.Cell0DsCoordinates(2, 2) = r / norm;
	mesh.Cell0DsCoordinates(0, 3) = r / norm;  mesh.Cell0DsCoordinates(1, 3) = -r / norm;  mesh.Cell0DsCoordinates(2, 3) = -r / norm;
	mesh.Cell0DsCoordinates(0, 4) = -r / norm;  mesh.Cell0DsCoordinates(1, 4) = r / norm;  mesh.Cell0DsCoordinates(2, 4) = r / norm;
	mesh.Cell0DsCoordinates(0, 5) = -r / norm;  mesh.Cell0DsCoordinates(1, 5) = r / norm;  mesh.Cell0DsCoordinates(2, 5) = -r / norm;
	mesh.Cell0DsCoordinates(0, 6) = -r / norm;  mesh.Cell0DsCoordinates(1, 6) = -r / norm;  mesh.Cell0DsCoordinates(2, 6) = r / norm;
	mesh.Cell0DsCoordinates(0, 7) = -r / norm;  mesh.Cell0DsCoordinates(1, 7) = -r / norm;  mesh.Cell0DsCoordinates(2, 7) = -r / norm;
	
	double norm1= sqrt(3.0+sqrt(5.0));
    mesh.Cell0DsCoordinates(0, 8) = 0;  mesh.Cell0DsCoordinates(1, 8) = (r * phi) / norm1;  mesh.Cell0DsCoordinates(2, 8) = phi / norm1;
	mesh.Cell0DsCoordinates(0, 9) = 0;  mesh.Cell0DsCoordinates(1, 9) = (r * phi) / norm1;  mesh.Cell0DsCoordinates(2, 9) = -phi / norm1;
	mesh.Cell0DsCoordinates(0, 10) = 0;  mesh.Cell0DsCoordinates(1, 10) = (-r * phi) / norm1;  mesh.Cell0DsCoordinates(2, 10) = phi / norm1;
	mesh.Cell0DsCoordinates(0, 11) = 0;  mesh.Cell0DsCoordinates(1, 11) = (-r * phi) / norm1;  mesh.Cell0DsCoordinates(2, 11) = -phi / norm1;
	mesh.Cell0DsCoordinates(0, 12) = phi / norm1;  mesh.Cell0DsCoordinates(1, 12) = 0;  mesh.Cell0DsCoordinates(2, 12) = (r * phi) / norm1;
	mesh.Cell0DsCoordinates(0, 13) = -phi / norm1;  mesh.Cell0DsCoordinates(1, 13) = 0;  mesh.Cell0DsCoordinates(2, 13) = (r * phi) / norm1;
	mesh.Cell0DsCoordinates(0, 14) = phi / norm1;  mesh.Cell0DsCoordinates(1, 14) = 0;  mesh.Cell0DsCoordinates(2, 14) = (-r * phi) / norm1;
	mesh.Cell0DsCoordinates(0, 15) = -phi / norm1;  mesh.Cell0DsCoordinates(1, 15) = 0;  mesh.Cell0DsCoordinates(2, 15) = (-r * phi) / norm1;
	mesh.Cell0DsCoordinates(0, 16) = (r * phi) / norm1;  mesh.Cell0DsCoordinates(1, 16) = phi / norm1;  mesh.Cell0DsCoordinates(2, 16) = 0;
	mesh.Cell0DsCoordinates(0, 17) = (r * phi) / norm1;  mesh.Cell0DsCoordinates(1, 17) = -phi / norm1;  mesh.Cell0DsCoordinates(2, 17) = 0;
	mesh.Cell0DsCoordinates(0, 18) = (-r * phi) / norm1;  mesh.Cell0DsCoordinates(1, 18) = phi / norm1;  mesh.Cell0DsCoordinates(2, 18) = 0;
	mesh.Cell0DsCoordinates(0, 19) = (-r * phi) / norm1;  mesh.Cell0DsCoordinates(1, 19) = -phi / norm1;  mesh.Cell0DsCoordinates(2, 19) = 0;

    mesh.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    
    // Lati
    mesh.Cell1DsId.reserve(30);
    Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

    mesh.Cell1DsExtrema = MatrixXi::Zero(30, 2);
    mesh.Cell1DsExtrema <<
        0, 1,  // Lato 0
        1, 2,  // Lato 1
        2, 3,  // Lato 2
        3, 4,  // Lato 3
        4, 0,  // Lato 4
        1, 5,  // Lato 5
        6, 5,  // Lato 6
        6, 7,  // Lato 7
        7, 2,  // Lato 8
        6, 10,  // Lato 9
        8, 9,  // Lato 10
        9, 10,  // Lato 11
        10, 11,  // Lato 12
        5, 12,  // Lato 13
        11, 12,  // Lato 14
        12, 13,  // Lato 15
        15, 11,  // Lato 16
        14, 15,  // Lato 17
        13, 14,  // Lato 18
        0, 13,  // Lato 19
        15, 16,  // Lato 20
        14, 18,  // Lato 21
        18, 4,  // Lato 22
        18, 17,  // Lato 23
        17, 16,  // Lato 24
        17, 19,  // Lato 25
        16, 9,  // Lato 26
        19, 8,  // Lato 27
        19, 3,  // Lato 28
        8, 7;  // Lato 29

    // Facce
    mesh.Cell2DsId.reserve(12);
    mesh.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    
    mesh.Cell2DsVertices.resize(12);
	for (auto& v : mesh.Cell2DsVertices)
		v.reserve(5);
		
	mesh.Cell2DsEdges.resize(12);
	for (auto& e : mesh.Cell2DsEdges)
		e.reserve(5);

     mesh.Cell2DsVertices = {
        {0, 1, 2, 3, 4},  // Faccia 0
        {1, 5, 6, 7, 2},  // Faccia 1
        {2, 7, 8, 19, 3},  // Faccia 2
        {10, 6, 7, 8, 9},  // Faccia 3
        {9, 8, 19, 17, 16},  // Faccia 4
        {3, 19, 17, 18, 4},  // Faccia 5
        {16, 17, 18, 14, 15},  // Faccia 6
        {15, 14, 13, 12, 11},   // Faccia 7
        {0, 4, 18, 14, 13},  // Faccia 8
        {13, 12, 5, 1, 0},  // Faccia 9
        {11, 12, 5, 6, 10},  // Faccia 10
        {9, 16, 15, 11, 10},  // Faccia 11
    };

    mesh.Cell2DsEdges = {
        {0, 1, 2, 3, 4},  // Faccia 0
        {5, 6, 7, 8, 1},  // Faccia 1
        {8, 29, 27, 28, 2},  // Faccia 2
        {9, 7, 29, 10, 11},  // Faccia 3
        {10, 27, 25, 24, 26},  // Faccia 4
        {28, 25, 23, 22, 3},  // Faccia 5
        {24, 23, 21, 17, 20},  // Faccia 6
        {17, 18, 15, 14, 16},   // Faccia 7
        {4, 22, 21, 18, 19},  // Faccia 8
        {15, 13, 5, 0, 19},  // Faccia 9
        {14, 13, 6, 9, 12},  // Faccia 10
        {26, 20, 16, 12, 11},  // Faccia 11
    };

    // Poliedro
    mesh.Cell3DsId.reserve(1);
    mesh.Cell3DsVertices.reserve(20);
    mesh.Cell3DsEdges.reserve(30);
    mesh.Cell3DsFaces.reserve(12);
    
    mesh.Cell3DsId = {0};
    mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5,6 ,7 ,8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    mesh.Cell3DsEdges = {0, 1, 2, 3,  4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
    mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    ExportDual(const PolyhedralMesh& mesh, PolyhedralMesh& dualMesh);
    CreateTxtFiles(const PolyhedralMesh& dualMesh);
    
    return true;
}

bool ExportIcosahedron(PolyhedralMesh& mesh) {

    // Vertici
    double r = 1.0;
    double phi = (1.0 + sqrt(5.0)) / 2.0;

    mesh.Cell0DsCoordinates = MatrixXd::Zero(12, 3);
    mesh.Cell0DsId.reserve(12);
	
	double norm = sqrt(10.0+2.0*sqrt(5.0))/2.0;
    mesh.Cell0DsCoordinates(0, 0) = 0.0; mesh.Cell0DsCoordinates(1, 0) = r / norm; mesh.Cell0DsCoordinates(2, 0) = phi / norm;
    mesh.Cell0DsCoordinates(0, 1) = 0.0; mesh.Cell0DsCoordinates(1, 1) = -r / norm; mesh.Cell0DsCoordinates(2, 1) = phi / norm;
    mesh.Cell0DsCoordinates(0, 2) = 0.0; mesh.Cell0DsCoordinates(1, 2) = r / norm; mesh.Cell0DsCoordinates(2, 2) = -phi / norm;
    mesh.Cell0DsCoordinates(0, 3) = 0.0; mesh.Cell0DsCoordinates(1, 3) = -r / norm; mesh.Cell0DsCoordinates(2, 3) = -phi / norm;
    mesh.Cell0DsCoordinates(0, 4) = r / norm; mesh.Cell0DsCoordinates(1, 4) = phi/ norm; mesh.Cell0DsCoordinates(2, 4) = 0.0;
    mesh.Cell0DsCoordinates(0, 5) = -r / norm; mesh.Cell0DsCoordinates(1, 5) = phi / norm; mesh.Cell0DsCoordinates(2, 5) = 0.0;
    mesh.Cell0DsCoordinates(0, 6) = r / norm; mesh.Cell0DsCoordinates(1, 6) = -phi / norm; mesh.Cell0DsCoordinates(2, 6) = 0.0;
    mesh.Cell0DsCoordinates(0, 7) = -r / norm; mesh.Cell0DsCoordinates(1, 7) = -phi / norm; mesh.Cell0DsCoordinates(2, 7) = 0.0;
    mesh.Cell0DsCoordinates(0, 8) = phi / norm; mesh.Cell0DsCoordinates(1, 8) = 0.0; mesh.Cell0DsCoordinates(2, 8) = r / norm;
    mesh.Cell0DsCoordinates(0, 9) = -phi / norm; mesh.Cell0DsCoordinates(1, 9) = 0.0; mesh.Cell0DsCoordinates(2, 9) = r / norm;
    mesh.Cell0DsCoordinates(0, 10) = phi / norm; mesh.Cell0DsCoordinates(1, 10) = 0.0; mesh.Cell0DsCoordinates(2, 10) = -r / norm;
    mesh.Cell0DsCoordinates(0, 11) = -phi / norm; mesh.Cell0DsCoordinates(1, 11) = 0.0; mesh.Cell0DsCoordinates(2, 11) = -r / norm;
	
	
    mesh.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    // Lati
    mesh.Cell1DsId.reserve(30);
    mesh.Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

    mesh.Cell1DsExtrema = matrixXi::Zero(30, 2);
    mesh.Cell1DsExtrema <<
    	0, 1,  // Lato 0
        1, 2,  // Lato 1
        2, 0,  // Lato 2
        3, 0,  // Lato 3
        4, 0,  // Lato 4
        0, 5,  // Lato 5
        1, 5,  // Lato 6
        3, 2,  // Lato 7
        5, 4,  // Lato 8
        4, 3,  // Lato 9
        4, 7,  // Lato 10
        5, 7,  // Lato 11
        5, 6,  // Lato 12
        1, 6,  // Lato 13
        1, 10,  // Lato 14
        4, 8,  // Lato 15
        2, 10,  // Lato 16
        3, 8,  // Lato 17
        2, 9,  // Lato 18
        6, 10,  // Lato 19
        10, 9,  // Lato 20
        6, 7,  // Lato 21
        7, 8,  // Lato 22
        8, 9,  // Lato 23
        7, 11,  // Lato 24
        8, 11,  // Lato 25
        10, 11,  // Lato 26
        6, 11,  // Lato 27
        9, 11,  // Lato 28
        3, 9;  // Lato 29

	// Facce
	mesh.Cell2DsId.reserve(20);
	mesh.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
	
	mesh.Cell2DsVertices.resize(20);
	for (auto& v : mesh.Cell2DsVertices)
		v.reserve(3);
		
	mesh.Cell2DsEdges.resize(20);
	for (auto& e : mesh.Cell2DsEdges)
		e.reserve(3);

	mesh.Cell2DsVertices = {
		{0, 1, 2},  // Faccia 0
        {0, 2, 3},  // Faccia 1
        {0, 3, 4},  // Faccia 2
        {0, 4, 5},  // Faccia 3
        {0, 5, 1},  // Faccia 4
        {5, 7, 4},  // Faccia 5
        {4, 7, 8},  // Faccia 6
        {4, 8, 3},   // Faccia 7
        {3, 8, 9},  // Faccia 8
        {9, 3, 2},  // Faccia 9
        {2, 10, 9},  // Faccia 10
        {2, 1, 10},  // Faccia 11
        {1, 6, 10},  // Faccia 12
        {1, 5, 6},  // Faccia 13
        {5, 6, 7},  // Faccia 14
        {6, 7, 11},  // Faccia 15
        {7, 11, 8},  // Faccia 16
        {8, 11, 9},  // Faccia 17
        {9, 10, 11},  // Faccia 18
        {6, 10, 11},   // Faccia 19
    };
    mesh.Cell2DsEdges.reserve(20,3);
	mesh.Cell2DsEdges = {
		{0, 1, 2},  // Faccia 0
        {2, 7, 3},  // Faccia 1
        {3, 9, 4},  // Faccia 2
        {4, 8, 5},  // Faccia 3
        {5, 6, 0},  // Faccia 4
        {11, 10, 8},  // Faccia 5
        {10, 22, 15},  // Faccia 6
        {15, 17, 9},   // Faccia 7
        {17, 23, 29},  // Faccia 8
        {29, 7, 18},  // Faccia 9
        {16, 20, 18},  // Faccia 10
        {1, 14, 16},  // Faccia 11
        {13, 19, 14},  // Faccia 12
        {6, 12, 13},  // Faccia 13
        {12, 21, 11},  // Faccia 14
        {21, 24, 27},  // Faccia 15
        {24, 25, 22},  // Faccia 16
        {25, 28, 23},  // Faccia 17
        {20, 26, 28},  // Faccia 18
        {19, 26, 27},   // Faccia 19
    };
    
    // Poliedro
    mesh.Cell3DsId.reserve(1);
    mesh.Cell3DsVertices.reserve(12);
    mesh.Cell3DsEdges.reserve(30);
    mesh.Cell3DsFaces.reserve(20);
    
    mesh.Cell3DsId = {0};
    mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5,6 ,7 ,8, 9, 10, 11};
    mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
	    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
    mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

    CreateTxtFiles(const PolyhedralMesh& mesh);
      
    return true;
}
}