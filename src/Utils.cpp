#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include "PolyhedralTriangulation.hpp"
#include <Eigen/Dense>
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralTriangulation; 

namespace PolyhedralLibrary
{
Vector3i ComputeVEF(unsigned int q, int b, int c)
{
    Vector3i VEF;  // inizializzo un vettore nullo di lunghezza 3
	unsigned int V, E, F;
	
	if (b!=c){
		unsigned int T = b * b + b * c + c * c;

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
	}
	
	if (b==c){
		unsigned int numV, numE, numF;
		if (q==3){
			numV = 4;
			numE = 6;
			numF = 4;
		}
		if (q==4){
			numV = 6;
			numE = 12;
			numF = 8;
		}
		if (q==5){
			numV = 12;
			numE = 30;
			numF = 20;
		}
		V = numV + numE*(2*b-1)+numF*((3*b*b)/2 - (3*b)/2 +1);
		E = numE*(2*b)+numF*((9*b*b)/2 + (3*b)/2);
		F = numF*(3*b*b + 3*b);	
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
bool EdgeIsDupe(const PolyhedralMesh& mesh, const Vector2i& e){
    	Vector2i v=e;
        for(size_t i=0; i<mesh.Cell1DsId.size(); i++){
            if(mesh.Cell1DsExtrema.col(i)==v)
                return true;
            swap(v[0],v[1]);
            if(mesh.Cell1DsExtrema.col(i)==v)
                return true;
        }
        return false;
    }


void OrderFaces(const vector<int>& unordered_faces, vector<int>& ordered_faces, const PolyhedralMesh& Mesh) {	
			//Il vettore di facce rimanenti contiene le facce ancora da ordinare, esso è all'inizio uguale a tutto il vettore.
			vector<int> remaining_faces = unordered_faces;

			// Inizio da una faccia qualsiasi, e la metto nel vettore ordinato
			int current = remaining_faces[0];
			ordered_faces.push_back(current);
			remaining_faces.erase(remaining_faces.begin());

			//Fino a quando ho facce da ordinare
			while (!remaining_faces.empty()) {
				//Prendo gli edge della faccia corrente
				vector<unsigned int> current_edges = Mesh.Cell2DsEdges[current];
				bool found = false;
				int found_index;

				// Cerco tra le facce da ordinare quella che ha un edge in comune con la corrente
				for (size_t i = 0; i < remaining_faces.size(); i++) {
					vector<unsigned int> candidate_edges = Mesh.Cell2DsEdges[remaining_faces[i]];

					// Verifico la presenza dell' edge in comune
					for (unsigned int e1 : current_edges) {
						for (unsigned int e2 : candidate_edges) {
							if (e1 == e2) {
								found = true;
								found_index = i;
								break;
							}
						}
						if (found) 
							break;
					}
					if (found) 
						break;
				}

				if (found) {
					//Se ho trovato l'edge comune, allora metto la faccia con l'edge comune dentro a ordered_faces
					//e aggiorno current, in modo da ripartire con il while cercando la faccia che ha un edge in comune con
					//quella appena trovata
					current = remaining_faces[found_index];
					ordered_faces.push_back(current);
					remaining_faces.erase(remaining_faces.begin() + found_index);
				}
			}
		}

bool GenerateDual(const PolyhedralMesh& baseMesh, PolyhedralMesh& dualMesh) {
			int baricenter_id = 0;
			int face_id = 0;
			int edge_id = 0;
			
			//Il Duale ha un numero di facce uguale al numero di vertici del poliedro di partenza
			dualMesh.NumCell0Ds = baseMesh.NumCell2Ds;
			
			//Il Duale ha lo stesso numero di edges del poliedro di partenza, grazie alla formula di Eulero, qualunque sia la varietà su cui si fa la mesh :)
			dualMesh.NumCell1Ds = baseMesh.NumCell1Ds;
			
			//Il Duale ha un numero di facce uguale al numero di vertici del poliedro di partenza
			dualMesh.NumCell2Ds = baseMesh.NumCell0Ds;
			dualMesh.Cell0DsId.reserve(dualMesh.NumCell0Ds);
			dualMesh.Cell0DsCoordinates = MatrixXd::Zero(3,dualMesh.NumCell0Ds);	
			
			dualMesh.Cell1DsExtrema = MatrixXi::Zero(2, dualMesh.NumCell1Ds);
			dualMesh.Cell1DsId.reserve(dualMesh.NumCell1Ds);
			
			dualMesh.Cell2DsId.reserve(dualMesh.NumCell2Ds);
			dualMesh.Cell2DsEdges.resize(dualMesh.NumCell2Ds);
			dualMesh.Cell2DsNumEdges.resize(dualMesh.NumCell2Ds);
			dualMesh.Cell2DsVertices.resize(dualMesh.NumCell2Ds);
			dualMesh.Cell2DsNumVertices.resize(dualMesh.NumCell2Ds);
			int duplicate_id = 0;
			//Mappa che associa all'id della faccia l'id del baricentro corrispondente, da usare se cambiassimo gli id dopo, per ora sono uguali
			map<int, int> Faces_bar;
			for (const auto& id : baseMesh.Cell2DsId) {
				Vector3d baricenter;
				
				// salvo in 3 vettori le coordinate dei vertici della faccia corrente del solido platonico
				Vector3d Vertex1 = baseMesh.Cell0DsCoordinates.col(baseMesh.Cell2DsVertices[id][0]);
				Vector3d Vertex2 = baseMesh.Cell0DsCoordinates.col(baseMesh.Cell2DsVertices[id][1]);
				Vector3d Vertex3 = baseMesh.Cell0DsCoordinates.col(baseMesh.Cell2DsVertices[id][2]);
				
				// salvo le coordinate del baricentro
				baricenter = (1.0/3.0)*Vertex1 + (1.0/3.0)*Vertex2 + (1.0/3.0)*Vertex3;
				baricenter = baricenter/baricenter.norm();
				dualMesh.Cell0DsId.push_back(baricenter_id);
				// salvo le coordinate dentro alle Cell0DsCoordinates del poliedro Duale
				dualMesh.Cell0DsCoordinates(0, id) = baricenter(0);
				dualMesh.Cell0DsCoordinates(1,id) = baricenter(1);
				dualMesh.Cell0DsCoordinates(2, id) = baricenter(2);
				
				//Associo all'id della faccia l'id del baricentro nella mappa, per ora ridondante ma poi sarà meglio
				Faces_bar[id] = baricenter_id;
				baricenter_id ++;
			}

			for(const auto& vertex_id: baseMesh.Cell0DsId){
				
				//Vettore che contiene le facce che hanno il vertice in comune
				vector<int> VertexFaces;
				for(const auto& face_id: baseMesh.Cell2DsId){
					for(int j = 0; j < 3; j++){
						if (baseMesh.Cell2DsVertices[face_id][j] == vertex_id){
							//Se la faccia a cui sono arrivato contiene il vertice, allora la aggiungo al vettore 
							//delle facce avente quel vertice in comune
							VertexFaces.push_back(face_id);
							break;
						}
					}
				}
				
				//PROBLEMA: il vettore VertexFaces contiene tutte le facce adiacenti a un vertice ma NON è ordinato in modo sensato, 
				//per costruire gli edges devo ordinarlo in modo che ogni faccia nel vettore abbia come successiva la faccia adiacente, 
				//ovvero quella che ha l'edge in comune con la faccia corrente, per ordinare questo vettore chiamo la funzione OrderFaces
				//il vettore ordered_Faces è passato per riferimento, in modo che venga aggiornato dalla funzione order_Faces.
				vector<int> ordered_faces;
				OrderFaces(VertexFaces, ordered_faces, baseMesh);
				
				//la valenza del vertice è pari alla lunghezza del vettore di facce che condividono il vertice dato 
				//ATTENZIONE: Questa parte non è superflua, perché le valenze NON SONO sempre 3 per il generico solido geodetico!!!
				int valence = ordered_faces.size();
				vector<unsigned int> New_vertices;
				
				//qui associo a ogni faccia del poliedro di partenza l'id del vertice nel duale corrispondente
				for(const auto& VertexFace: ordered_faces)
					New_vertices.push_back(Faces_bar[VertexFace]);
				
				
				dualMesh.Cell2DsId.push_back(face_id);
				dualMesh.Cell2DsVertices[face_id] = New_vertices;
				dualMesh.Cell2DsEdges[face_id].resize(valence);
				
				dualMesh.Cell2DsNumVertices[face_id] = valence;
				dualMesh.Cell2DsNumEdges[face_id] = valence;
				
				//Questa è la creazione degli edges, praticamente identica a quella per il solido geodetico, 
				//con la differenza che il vettore di vertici della faccia ha tanti elementi quanti la valenza del 
				//vertice, e a parte questo dettaglio è tutto uguale!
				for (int k = 0; k < valence; k++) {
					int originVertex = dualMesh.Cell2DsVertices[face_id][k];
					int endVertex;
					if ( k == valence-1 )
						endVertex = dualMesh.Cell2DsVertices[face_id][0];
					else
						endVertex = dualMesh.Cell2DsVertices[face_id][k+1];
					Vector2i extrema(originVertex,endVertex);
					if(!EdgeIsDupe(dualMesh, extrema)){
						dualMesh.Cell1DsId.push_back(edge_id);
						dualMesh.Cell1DsExtrema(0, edge_id) = originVertex;
						dualMesh.Cell1DsExtrema(1, edge_id) = endVertex;
						dualMesh.Cell2DsEdges[face_id][k] = edge_id;
						edge_id++;
					}
					else
						dualMesh.Cell2DsEdges[face_id][k] = dualMesh.Cell1DsId[duplicate_id];
				}
				face_id++;
			}
			

			
			dualMesh.Cell1DsExtrema.conservativeResize(2, dualMesh.NumCell1Ds);
			dualMesh.Cell2DsNumVertices.resize(dualMesh.NumCell2Ds);
			dualMesh.Cell2DsNumEdges.resize(dualMesh.NumCell2Ds);
			dualMesh.Cell2DsVertices.resize(dualMesh.NumCell2Ds);
			dualMesh.Cell2DsEdges.resize(dualMesh.NumCell2Ds);
			
			// GENERAZIONE POLIEDRO
			dualMesh.NumCell3Ds++;
			dualMesh.Cell3DsId = {0};
			dualMesh.Cell3DsVertices = dualMesh.Cell0DsId;
			dualMesh.Cell3DsEdges = dualMesh.Cell1DsId;
			dualMesh.Cell3DsFaces = dualMesh.Cell2DsId;
	
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
	const double phi = (1.0 + sqrt(5.0)) / 2.0;
    const double r = 1.0/(sqrt(1 + pow(phi,2)));
	const double s = phi / (sqrt(1 + pow(phi,2)));

    mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 12);
    mesh.Cell0DsId.reserve(12);
	
    mesh.Cell0DsCoordinates(0, 0) = -r; mesh.Cell0DsCoordinates(1, 0) = s; mesh.Cell0DsCoordinates(2, 0) = 0.0;
    mesh.Cell0DsCoordinates(0, 1) = r; mesh.Cell0DsCoordinates(1, 1) = -s; mesh.Cell0DsCoordinates(2, 1) = 0.0;
    mesh.Cell0DsCoordinates(0, 2) = r; mesh.Cell0DsCoordinates(1, 2) = s; mesh.Cell0DsCoordinates(2, 2) = 0.0;
    mesh.Cell0DsCoordinates(0, 3) = -r; mesh.Cell0DsCoordinates(1, 3) = -s; mesh.Cell0DsCoordinates(2, 3) = 0.0;
    mesh.Cell0DsCoordinates(0, 4) = 0.0; mesh.Cell0DsCoordinates(1, 4) = r; mesh.Cell0DsCoordinates(2, 4) = s;
    mesh.Cell0DsCoordinates(0, 5) = 0.0; mesh.Cell0DsCoordinates(1, 5) = -r; mesh.Cell0DsCoordinates(2, 5) = -s;
    mesh.Cell0DsCoordinates(0, 6) = 0.0; mesh.Cell0DsCoordinates(1, 6) = -r; mesh.Cell0DsCoordinates(2, 6) = s;
    mesh.Cell0DsCoordinates(0, 7) = 0.0; mesh.Cell0DsCoordinates(1, 7) = r; mesh.Cell0DsCoordinates(2, 7) = -s;
    mesh.Cell0DsCoordinates(0, 8) = -s; mesh.Cell0DsCoordinates(1, 8) = 0.0; mesh.Cell0DsCoordinates(2, 8) = r;
    mesh.Cell0DsCoordinates(0, 9) = -s; mesh.Cell0DsCoordinates(1, 9) = 0.0; mesh.Cell0DsCoordinates(2, 9) = -r;
    mesh.Cell0DsCoordinates(0, 10) = s; mesh.Cell0DsCoordinates(1, 10) = 0.0; mesh.Cell0DsCoordinates(2, 10) = r;
    mesh.Cell0DsCoordinates(0, 11) = s; mesh.Cell0DsCoordinates(1, 11) = 0.0; mesh.Cell0DsCoordinates(2, 11) = -r;
	
	
    mesh.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

    // Lati
    mesh.Cell1DsId.reserve(30);
    mesh.Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

    mesh.Cell1DsExtrema = MatrixXi::Zero(2, 30);
    mesh.Cell1DsExtrema <<
    	0, 2,  // Lato 0
        0, 4,  // Lato 1
        0, 7,  // Lato 2
        0, 8,  // Lato 3
        0, 9,  // Lato 4
        1, 3,  // Lato 5
        1, 6,  // Lato 6
        1, 5,  // Lato 7
        1, 10,  // Lato 8
        1, 11,  // Lato 9
        2, 4,  // Lato 10
        2, 7,  // Lato 11
        2, 11,  // Lato 12
        2, 10,  // Lato 13
        3, 5,  // Lato 14
        3, 9,  // Lato 15
        3, 8,  // Lato 16
        3, 6,  // Lato 17
        4, 10,  // Lato 18
        4, 6,  // Lato 19
        4, 8,  // Lato 20
        5, 11,  // Lato 21
        5, 7,  // Lato 22
        5, 9,  // Lato 23
        6, 8,  // Lato 24
        6, 10,  // Lato 25
        7, 9,  // Lato 26
        7, 11,  // Lato 27
        8, 9,  // Lato 28
        10, 11;  // Lato 29

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
		{0, 2, 4},  // Faccia 0
        {0, 4, 8},  // Faccia 1
        {0, 8, 9},  // Faccia 2
        {0, 7, 9},  // Faccia 3
        {0, 2, 7},  // Faccia 4
        {2, 7, 11},  // Faccia 5
        {5, 7, 11},  // Faccia 6
        {5, 7, 9},   // Faccia 7
        {3, 5, 9},  // Faccia 8
        {3, 8, 9},  // Faccia 9
        {3, 6, 8},  // Faccia 10
        {4, 6, 8},  // Faccia 11
        {4, 6, 10},  // Faccia 12
        {1, 6, 10},  // Faccia 13
        {1, 3, 6},  // Faccia 14
        {1, 3, 5},  // Faccia 15
        {1, 5, 11},  // Faccia 16
        {1, 10, 11},  // Faccia 17
        {2, 10, 11},  // Faccia 18
        {2, 4, 10},   // Faccia 19
    };
    mesh.Cell2DsEdges.reserve(20);
	for (auto& edgeList : mesh.Cell2DsEdges) {
		edgeList.resize(3);
	}	
	mesh.Cell2DsEdges = {
		{0, 10, 1},  // Faccia 0
        {1, 20, 3},  // Faccia 1
        {3, 28, 4},  // Faccia 2
        {2, 26, 4},  // Faccia 3
        {0, 11, 2},  // Faccia 4
        {11, 27, 12},  // Faccia 5
        {22, 27, 21},  // Faccia 6
        {22, 26, 23},   // Faccia 7
        {14, 23, 15},  // Faccia 8
        {16, 28, 15},  // Faccia 9
        {17, 24, 16},  // Faccia 10
        {19, 24, 20},  // Faccia 11
        {19, 25, 18},  // Faccia 12
        {6, 25, 8},  // Faccia 13
        {5, 17, 6},  // Faccia 14
        {5, 14, 7},  // Faccia 15
        {7, 21, 9},  // Faccia 16
        {8, 29, 9},  // Faccia 17
        {13, 29, 12},  // Faccia 18
        {10, 18, 13},   // Faccia 19
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
    else {
    	GenerateTriangulatedMesh2(mesh,triMesh,b,c,VEF);
		cout << "Fatta triangolazione II" << endl; 
		}
      
    return true;
}
	
	/////////////////////////////////////////////////////////////////////////////////////////////
	
/**bool ShortestPath(PolyhedralMesh& mesh, unsigned int id_vertice_1, unsigned int id_vertice_2, unsigned int num_lati_iniziali) {
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
        unsigned int w = cammino[i + 1];

        
        for (unsigned int j = num_lati_iniziali; j < mesh.NumCell1Ds; j++) {
            unsigned int v1 = mesh.Cell1DsExtrema(0, j);
            unsigned int v2 = mesh.Cell1DsExtrema(1, j);
            if ((v1 == u && v2 == w) || (v1 == w && v2 == u)) {
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

void ExportParaView(PolyhedralMesh& mesh, bool cammino){
	
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
		{
		utilities.ExportPolygons("./Cell2Ds.inp",
                             mesh.Cell0DsCoordinates,
                             mesh.Cell2DsVertices,
                             {});
		}
	} /* else{
		
			
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
}
