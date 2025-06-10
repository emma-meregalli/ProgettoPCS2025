#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include "PolyhedralTriangulation.hpp"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace PolyhedralTriangulation; 

namespace PolyhedralLibrary
{
Vector3i ComputeVEF(unsigned int q, int b, int c)
{
    Vector3i VEF;  // inizializzo un vettore nullo di lunghezza 3

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

/////////////////////////////////////////////////////////////////////////////////////////////

void CreateTxtFiles(const PolyhedralMesh& mesh) {
    // Creazione Cell0Ds.txt
    ofstream Cell0Ds("Cell0Ds.txt");
    Cell0Ds << "ID;x;y;z\n";
    for (size_t i = 0; i < mesh.Cell0DsId.size(); i++) {
        Cell0Ds << mesh.Cell0DsId[i] << ";" << mesh.Cell0DsCoordinates(0,i) 
			<< ";" << mesh.Cell0DsCoordinates(1,i) << ";" << mesh.Cell0DsCoordinates(2,i) << "\n";
    }
    Cell0Ds.close();

    // Creazione Cell1Ds.txt 
    ofstream Cell1Ds("Cell1Ds.txt");
    Cell1Ds << "ID;Origin;End\n";
    for (size_t i = 0; i < mesh.Cell1DsId.size(); i++) {
        Cell1Ds << mesh.Cell1DsId[i] << ";" << mesh.Cell1DsExtrema(0, i) << ";"
             << mesh.Cell1DsExtrema(1, i) << "\n";
    }
    Cell1Ds.close();
	cout<<"letsgooo"<<endl;
    // Creazione Cell2Ds.txt 
    ofstream Cell2Ds("Cell2Ds.txt");
    Cell2Ds << "ID;NumVertices;NumEdges;Vertices;Edges\n";
    cout<<mesh.Cell2DsId.size()<<endl;
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
    for (size_t i = 0; i < mesh.Cell3DsId.size(); ++i) {
	    Cell3Ds << mesh.Cell3DsId[i];
	
	    // Vertici
	    for (unsigned int v = 0; v < mesh.Cell3DsVertices.size(); v++) {
	        Cell3Ds << ";" << mesh.Cell3DsVertices[v];
	    }
	
	    // Lati
	    for (unsigned int e = 0; e < mesh.Cell3DsEdges.size(); e++) {
	        Cell3Ds << ";" << mesh.Cell3DsEdges[e];
	    }
	
	    // Facce
	    for (unsigned int f = 0; f < mesh.Cell3DsFaces.size(); f++) {
	        Cell3Ds << ";" << mesh.Cell3DsFaces[f];
	    }
	
	    Cell3Ds << "\n";
	}
    Cell3Ds.close();
} 

/////////////////////////////////////////////////////////////////////////////////////////////

bool GenerateDual(const PolyhedralMesh& mesh, PolyhedralMesh& dualMesh){
	MatrixXd barycenters = {};
	barycenters.resize(mesh.NumCell2Ds,3);
    Vector3d tmp = Vector3d::Zero();
    for (unsigned int i=0; i<mesh.NumCell2Ds; i++){
        for (unsigned int v : mesh.Cell2DsVertices[i]){
            tmp += mesh.Cell0DsCoordinates.row(v);
        }
        tmp /= mesh.Cell2DsVertices[i].size(); 
        barycenters.row(i) = tmp;
    }

    // Riempie Cell0Ds del duale con i baricentri
    int numDualVertices = barycenters.size();
    dualMesh.Cell0DsCoordinates.resize(numDualVertices,3);
    for (int i = 0; i < numDualVertices; i++) {
        dualMesh.Cell0DsCoordinates.row(i) = barycenters.row(i);
        dualMesh.Cell0DsId.push_back(i);
    }

    // Mappa: vertice originale -> facce adiacenti che lo contengono
    map<int, vector<int>> vertexToFaces;
    for (size_t f = 0; f < mesh.Cell2DsVertices.size(); f++) {
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
            for (auto it = rest.begin(); it != rest.end(); it++) {
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
        std::vector<unsigned int> orderedFacesUnsigned(orderedFaces.begin(), orderedFaces.end());
		dualMesh.Cell2DsVertices.push_back(orderedFacesUnsigned);
        dualMesh.Cell2DsId.push_back(faceId);

        // Costruzione degli spigoli per la faccia
        vector<unsigned int> faceEdges;
        for (size_t i = 0; i < orderedFaces.size(); i++) {
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
        faceId++;
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
    
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////

bool ExportTetrahedron(PolyhedralMesh& mesh, PolyhedralMesh& triMesh, const int& b, const int& c) {
	
	// Vertici
    double r = sqrt(3.0) / 3.0;
    
    mesh.Cell0DsCoordinates.resize(3,4);
    mesh.Cell0DsId.reserve(4);
	
    mesh.Cell0DsCoordinates(0, 0) = r;   mesh.Cell0DsCoordinates(1, 0) = r;   mesh.Cell0DsCoordinates(2, 0) = r;  
    mesh.Cell0DsCoordinates(0, 1) = -r;  mesh.Cell0DsCoordinates(1, 1) = -r;  mesh.Cell0DsCoordinates(2, 1) = r; 
    mesh.Cell0DsCoordinates(0, 2) = -r;  mesh.Cell0DsCoordinates(1, 2) = r;   mesh.Cell0DsCoordinates(2, 2) = -r;
    mesh.Cell0DsCoordinates(0, 3) = r;   mesh.Cell0DsCoordinates(1, 3) = -r;  mesh.Cell0DsCoordinates(2, 3) = -r;
    
    mesh.Cell0DsId = {0,1,2,3};
    
    // Lati
    mesh.Cell1DsId.reserve(6);
    mesh.Cell1DsId = {0,1,2,3,4,5};
	
    mesh.Cell1DsExtrema.resize(2, 6);
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
    
    Vector3i VEF = ComputeVEF(3,b,c);
    if(b!=c)
		GenerateTriangulatedMesh1(mesh,triMesh,b,c,VEF);
    else
    	GenerateTriangulatedMesh2(mesh,triMesh,b,c,VEF);
	
    return true;
}
	
/////////////////////////////////////////////////////////////////////////////////////////////
	
bool ExportOctahedron(PolyhedralMesh& mesh, PolyhedralMesh& triMesh, const int& b, const int& c) {
	
    // Vertici
    double r = 1.0; 

    mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 6);
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
    mesh.Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
	
    mesh.Cell1DsExtrema = MatrixXi::Zero(2, 12);
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
    
    Vector3i VEF = ComputeVEF(4,b,c);
    if(b!=c)
		GenerateTriangulatedMesh1(mesh,triMesh,b,c,VEF);
    else
    	GenerateTriangulatedMesh2(mesh,triMesh,b,c,VEF);
		
    return true;
}
	
/////////////////////////////////////////////////////////////////////////////////////////////

bool ExportIcosahedron(PolyhedralMesh& mesh, PolyhedralMesh& triMesh, const int& b, const int& c) {

    // Vertici
    double r = 1.0;
    double phi = (1.0 + sqrt(5.0)) / 2.0;

    mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 12);
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

    mesh.Cell1DsExtrema = MatrixXi::Zero(2, 30);
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
    mesh.Cell2DsEdges.reserve(20);
	for (auto& edgeList : mesh.Cell2DsEdges) {
		edgeList.resize(3);
	}	
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
    mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
    mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
    
    Vector3i VEF = ComputeVEF(5,b,c);
    if(b!=c)
		GenerateTriangulatedMesh1(mesh,triMesh,b,c,VEF);
    else
    	GenerateTriangulatedMesh2(mesh,triMesh,b,c,VEF);
      
    return true;
}
	
	/////////////////////////////////////////////////////////////////////////////////////////////
	
/*bool ShortestPath(PolyhedralMesh& mesh, unsigned int v1, unsigned int v2, unsigned int num_lati_iniziali) {
    if (id_vertice_1 >= mesh.NumCell0Ds || id_vertice_2 >= mesh.NumCell0Ds) {
        cerr << "ID dei vertici non valido." << endl;
        return false;
    }

    unsigned int N = mesh.NumCell0Ds;

    //Inizializzo la lista di adiacenza
	//Uso un vettore di vettori. L'i-esimo vettore corrisponde all'i-esimo vertice, ed è costituito da coppie vertice-peso, cioè contiene i vertici vicini e la relativa distanza
    vector<vector<pair<unsigned int, double>>> LA(N); 

	//Riempio la lista di adiacenza
    for (unsigned int i = num_lati_iniziali; i < mesh.NumCell1Ds; i++) {//Considero solo i nuovi lati, quelli prima della triangolazione "non esistono più"
	
        unsigned int v1 = mesh.Cell1DsExtrema(0, i);
        unsigned int v2 = mesh.Cell1DsExtrema(1, i);

        Eigen::Vector3d c1 = mesh.Cell0DsCoordinates.col(v1);
        Eigen::Vector3d c2 = mesh.Cell0DsCoordinates.col(v2);
        double peso = (c1 - c2).norm();

        LA[v1].push_back({v2, peso});
        LA[v2].push_back({v1, peso});
    }

    // Inizializzazione Dijkstra
    vector<double> distanze(N, numeric_limits<double>::infinity()); //Vettore delle distanze, inizializzate tutte ad un numero molto grande
    vector<int> pred(N, -1); //Vettore dei predecessori, inizializzati tutti a -1
    vector<bool> visitati(N, false); //Vettore dei vertici visitati, inizializzati a false

    distanze[id_vertice_1] = 0.0; //La distanza del primo vertice da sè stesso è nulla
    
    priority_queue<pair<double,unsigned int>, vector<pair<double, unsigned int>>, greater<pair<double, unsigned int>>> pq; //Inizializzo la coda di priorità
	//La pq è costituita da una coppia (pair) dove salviamo la distanza e l'identificativo del vertice. Il greater ci permette di ordinare la coda in modo da prendere sempre 
	//per primo il vertice con la distanza minore
	//All'interno della pq c'è un altro contenitore, ovvero un vettore anch'esso costituito da coppie distanza, vertice, e che ci permette di memorizzare gli elementi
	
    pq.push({0.0, id_vertice_1}); //Aggiungo il primo vertice
	
	
	//Uso Dijkstra per il cammino minimo pesato
	
    while (!pq.empty()) {
		//Estraggo l'elemento dalla coda
        unsigned int u = pq.top().second; 
        pq.pop();


        if (visitati[u]) continue;
        visitati[u] = true;

        for (auto& [v, peso] : LA[u]) { //Vado a prendere ogni vertice con il relativo peso nella lista di adiacenza
            if (visitati[v]) continue; //Se l'ho già visitato non lo considero più

            if (distanze[u] + peso < distanze[v]) { //Altrimenti, controllo se aggiungerlo alla coda
                distanze[v] = distanze[u] + peso;
                pred[v] = u;
                pq.push({distanze[v], v});
            }
        }
    }

    // Ricostruisco il cammino (controllo se esiste)
    if (pred[id_vertice_2] == -1) {
        cerr << "Nessun cammino trovato da " << id_vertice_1 << " a " << id_vertice_2 << endl;
        return false;
    }


    vector<unsigned int> cammino;
	int nodo_corrente = id_vertice_2;

	//Costruisco il cammino partendo dall'ultimo vertice
	//Faccio così perchè mi serve la condizione sul predecessore (pari a -1), in modo tale da fermarmi una volta arrivato al nodo iniziale
	while (nodo_corrente != -1) {
		cammino.push_back(nodo_corrente);
		nodo_corrente = pred[nodo_corrente];  
	}

    reverse(cammino.begin(), cammino.end()); //Inverto l'ordine in modo da avere il cammino vero

    //Inizializzo le proprietà della mesh
    mesh.Cell0DsShortPath.resize(mesh.NumCell0Ds, 0);
    mesh.Cell1DsShortPath.resize(mesh.NumCell1Ds, 0);
	
	
	//Aggiungo i vertici visitati
    for (int v : cammino) {
        mesh.Cell0DsShortPath[v] = 1;
    }

 
    double Lunghezza_tot = 0.0;
    unsigned int Lunghezza_cammino = 0;
	
	//Aggiungo i lati visitati, aggiorno il numero di lati e la distanza totale percorsa
    for (unsigned int i = 0; i < cammino.size() - 1; i++) {
        unsigned int u = cammino[i];
        unsigned int v = cammino[i + 1];

        
        for (unsigned int j = num_lati_iniziali; j < mesh.NumCell1Ds; j++) {
            unsigned int v1 = mesh.Cell1DsExtrema(0, j);
            unsigned int v2 = mesh.Cell1DsExtrema(1, j);
            if ((v1 == u && v2 == v) || (v1 == v && v2 == u)) {
                mesh.Cell1DsShortPath[j] = 1;
                Eigen::Vector3d c1 = mesh.Cell0DsCoordinates.col(v1);
                Eigen::Vector3d c2 = mesh.Cell0DsCoordinates.col(v2);
                Lunghezza_tot += (c1 - c2).norm();
                Lunghezza_cammino++;
                break;
            }
        }
		
    }
	


    cout << "Numero di archi nel cammino: " << Lunghezza_cammino << endl;
    cout << "Lunghezza totale: " << Lunghezza_tot << endl;
	
	return true;
}*/

////////////////////////////////////////////////////////////////////////

/*void ExportParaView(PolygonalMesh& mesh, bool cammino){
	
	if(!cammino){
		
		Gedim::UCDUtilities utilities;
		{ 
		utilities.ExportPoints("./Cell0Ds.inp",
                               mesh.Cell0DsCoordinates);
		}
		{
		utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema,
                                 {});
		}
	}else{
		
			
		Gedim::UCDUtilities utilities;
		{
		vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

		cell0Ds_properties[0].Label = "ShortPathVertici";
		cell0Ds_properties[0].UnitLabel = "-";
		cell0Ds_properties[0].NumComponents = 1;

		vector<double> ShortPath_v(mesh.NumCell0Ds, 0.0);
		for (unsigned int i = 0; i < mesh.NumCell0Ds; ++i) {
			ShortPath_v[i] = static_cast<double>(mesh.Cell0DsShortPath[i]);
		}

		cell0Ds_properties[0].Data = ShortPath_v.data();

		utilities.ExportPoints("./Cell0Ds.inp",
							   mesh.Cell0DsCoordinates,
							   cell0Ds_properties);
		}
			
		{
		std::vector<Gedim::UCDProperty<double>> cell1Ds_properties(2);

		cell1Ds_properties[0].Label = "ShortPathLati";
		cell1Ds_properties[0].UnitLabel = "-";
		cell1Ds_properties[0].NumComponents = 1;

		std::vector<double> ShortPath_e(mesh.NumCell1Ds);
		for (unsigned int i = 0; i < mesh.NumCell1Ds; ++i)
			ShortPath_e[i] = static_cast<double>(mesh.Cell1DsShortPath[i]);
		cell1Ds_properties[0].Data = ShortPath_e.data();
		
		

		cell1Ds_properties[1].Label = "Lati Esistenti";
		cell1Ds_properties[1].UnitLabel = "-";
		cell1Ds_properties[1].NumComponents = 1;

		std::vector<double> LatiEsistenti(mesh.NumCell1Ds);
		for (unsigned int i = 0; i < mesh.NumCell1Ds; ++i)
			LatiEsistenti[i] = static_cast<double>(mesh.Cell1DsEsistente[i]);
		cell1Ds_properties[1].Data = LatiEsistenti.data();



		utilities.ExportSegments("./Cell1Ds.inp",
								mesh.Cell0DsCoordinates,
								mesh.Cell1DsExtrema,
								{},
								cell1Ds_properties);
		}
	}
}*/

}
