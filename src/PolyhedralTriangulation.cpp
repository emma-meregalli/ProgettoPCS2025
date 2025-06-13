#include <vector>  
#include <array>  
#include <Eigen/Dense>  
#include <limits> 
#include <utility>       
#include "Utils.hpp"   
#include "PolyhedralMesh.hpp" 

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

namespace PolyhedralTriangulation {

    bool VertexIsDupe(const PolyhedralMesh& mesh, const Vector3d& v, unsigned int& original_id){
        //Fisso una tolleranza per confrontare i vertici
	    double tol=1e-12;
        // Confronto con tutti i vertici già inseriti nella lista
        for (size_t i = 0; i < mesh.Cell0DsId.size(); i++) {
            if ((mesh.Cell0DsCoordinates.col(i) - v).norm() < tol){    //Se il vertice esiste, allora restituisco il suo ID (bisogna fare un controllo con la tolleranza?)
				original_id=i;
				return true;  
            }
        }
        return false;
    }

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
    
    // Funzione principale per triangolare le facce e popolare la mesh triangolata
    bool GenerateTriangulatedMesh1(
        PolyhedralMesh& baseMesh,     // Mesh di partenza (con facce non triangolate)
        PolyhedralMesh& triMesh,      // Mesh di output triangolata
        const unsigned int& b, const unsigned int& c, // Parametri della suddivisione
        const Vector3i& triDimensions) // Dimensione di (V,E,F) della mesh triangolata (con duplicati)
    {
        unsigned int level = b + c; // Numero di suddivisioni laterali per triangolo

        // Inizializzazione della struttura dati della mesh triangolata

        // Allocazione memoria per vertici (0D)
		triMesh.Cell0DsCoordinates = MatrixXd::Zero(3,triDimensions[0]); 
		triMesh.Cell0DsId.reserve(triDimensions[0]);
		
        // Allocazione memoria per lati (1D)
        triMesh.Cell1DsId.reserve(triDimensions[1]);
        triMesh.Cell1DsExtrema = MatrixXi::Zero(2, triDimensions[1]);
	
        // Allocazione memoria per facce (2D)
        triMesh.Cell2DsId.reserve(triDimensions[2]);
        triMesh.Cell2DsEdges.reserve(triDimensions[2]);
        triMesh.Cell2DsVertices.reserve(triDimensions[2]);
        for (auto& edgeList : triMesh.Cell2DsEdges) {
            edgeList.resize(3); // Ogni faccia ha 3 spigoli
        }
		for (auto& vertList : triMesh.Cell2DsVertices) {
        vertList.resize(3); // Ogni faccia ha 3 vertici
		}	
		
        // Id per vertici, lati e facce
        unsigned int vCount = 0;
		unsigned int eCount = 0;
		unsigned int fCount = 0;
		
		// faceIdx = 1, i = 1, j = 0

        // Ciclo su tutte le facce della mesh di base
        for (unsigned int faceIdx = 0; faceIdx < baseMesh.Cell2DsId.size(); faceIdx++) 
        {
            const auto& faceVerts = baseMesh.Cell2DsVertices[faceIdx]; //Prendo i tre vertici della faccia corrente
        
            // Coordinate dei 3 vertici del triangolo originale
            Vector3d A = baseMesh.Cell0DsCoordinates.col(faceVerts[0]); // Vertice A
            Vector3d B = baseMesh.Cell0DsCoordinates.col(faceVerts[1]); // Vertice B
            Vector3d C = baseMesh.Cell0DsCoordinates.col(faceVerts[2]); // Vertice C
			
            vector<vector<unsigned int>> grid; // Griglia di vertici interni alla faccia
            // la griglia ha level = b + c righe e ogni riga i ha i + 1 elementi (forma triangolare)
            // Costruzione della griglia interplata sulla faccia
            for (unsigned int i = 0; i <= level; i++){
                vector<unsigned int> row; //riga corrente dei vertici 
				unsigned int original_id;
                // Calcolo il punto iniziale e finale della riga i-esima 
                // Partiziono il lato in base al valore di b e c
                Vector3d from = ((double)i / level) * B + ((double)(level - i) / level) * A;
                Vector3d to = ((double)i / level) * C + ((double)(level - i) / level) * A;
                for (unsigned int j = 0; j <= i; j++) { 
                    // Interpolo tra from e to per ottenere un punto interno
                    Vector3d pos;
                    if (i == 0) {
                        pos = A;
                    } else {
                        pos = ((double)j / i) * to + ((double)(i - j) / i) * from;
                    }
                    pos = pos/pos.norm();
					if(!VertexIsDupe(triMesh, pos, original_id)){
						triMesh.Cell0DsId.push_back(vCount);           // Salva ID
						for(unsigned int n=0; n<3; n++){
							triMesh.Cell0DsCoordinates(n,vCount) = pos(n); // Salva posizione
						}
                    	row.push_back(vCount); // Aggiungi indice del vertice alla riga corrente
                    	vCount++; // Avanza contatore vertice
					}
					else{
						row.push_back(original_id);
					}
                }
                grid.push_back(row); // Aggiungi riga alla griglia
            }

            //Creiamo i nuovi lati dati dalla triangolazione e aggiorniamo la lista dei lati
            for(size_t i=0; i<grid.size(); i++){
                Vector2i extrema;
                for(size_t j=0; j<=i; j++){
                    if(i<grid.size()-1){
                    	extrema[0] = grid[i][j];
						extrema[1] = grid[i + 1][j];
                    	
                    	if(!EdgeIsDupe(triMesh, extrema)){
                    		triMesh.Cell1DsId.push_back(eCount);
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  //lato sotto a sinistra
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i + 1][j];
	                        eCount++;	
						}
						
						extrema[0] = grid[i][j];
						extrema[1] = grid[i + 1][j+1];
						if(!EdgeIsDupe(triMesh, extrema)){
	                        triMesh.Cell1DsId.push_back(eCount);
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  //lato sotto a destra
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i + 1][j + 1];
	                        eCount++;
	                    }
                    }
                    
                    if(i>0 && j<=i-1){
                    	extrema[0] = grid[i][j];
						extrema[1] = grid[i][j+1];
                    	if(!EdgeIsDupe(triMesh, extrema)){
	                        triMesh.Cell1DsId.push_back(eCount);
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  //lato accanto a destra
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i][j + 1];
	                        eCount++;
                    	}
                    }
                }
            }

			//Creiamo le nuove facce dopo la triangolazione
			for (size_t i=0; i<grid.size()-1; i++){
				for (size_t j=0; j<=i; j++){
					vector<unsigned int> v1 = {grid[i][j],grid[i+1][j],grid[i+1][j+1]};
					vector<unsigned int> e1;
					for(unsigned int v=0; v<3; v++){
						unsigned int from=v1[v];
						unsigned int to= v1[(v + 1) % 3];
						for(unsigned int k=0; k< triMesh.Cell1DsId.size(); k++){
							if((from==triMesh.Cell1DsExtrema(0, k) && to==triMesh.Cell1DsExtrema(1, k)) || (from==triMesh.Cell1DsExtrema(1, k) && to==triMesh.Cell1DsExtrema(0, k))){
								e1.push_back(k);
								break;
							}
						}
					}
					triMesh.Cell2DsId.push_back(fCount);
					triMesh.Cell2DsVertices.push_back(v1);
					triMesh.Cell2DsEdges.push_back(e1);
					fCount++;
					
					if(i>0 && j<=i-1){
						vector<unsigned int> v2 = {grid[i][j],grid[i][j+1],grid[i+1][j+1]};
						vector<unsigned int> e2;
						for(unsigned int v=0; v<3; v++){
							unsigned int from=v2[v];
							unsigned int to=v2[(v+1)%3];
							for(unsigned int k=0; k< triMesh.Cell1DsId.size(); k++){
								if((from==triMesh.Cell1DsExtrema(0, k) && to==triMesh.Cell1DsExtrema(1, k)) || (from==triMesh.Cell1DsExtrema(1, k) && to==triMesh.Cell1DsExtrema(0, k))){
									e2.push_back(k);
									break;
								}
							}	
						}
						triMesh.Cell2DsId.push_back(fCount);
						triMesh.Cell2DsVertices.push_back(v2);
						triMesh.Cell2DsEdges.push_back(e2);
						fCount++;
					}		
			    }
		    }
		}
    
		//Creiamo il poliedro triangolato
	    triMesh.Cell3DsId = {0};
		triMesh.Cell3DsVertices = triMesh.Cell0DsId;
	    triMesh.Cell3DsEdges = triMesh.Cell1DsId;
	    triMesh.Cell3DsFaces = triMesh.Cell2DsId;
		
		// Aggiorna i conteggi delle celle
	    triMesh.NumCell0Ds = triMesh.Cell0DsId.size();
	    triMesh.NumCell1Ds = triMesh.Cell1DsId.size();
		triMesh.NumCell2Ds = triMesh.Cell2DsId.size();
	    triMesh.NumCell3Ds = 1;

		return true;  
	}
	
	// Funzione per triangolare le facce e popolare la mesh triangolata di classe II
bool GenerateTriangulatedMesh2(
    PolyhedralMesh& baseMesh,
    PolyhedralMesh& triMesh,
    const unsigned int& b, const unsigned int& c, // Parametri della suddivisione
    const Vector3i& triDimensions) // Dimensione di (V,E,F) della mesh triangolata
	{
	    unsigned int level = b + c;
	    triMesh.Cell0DsCoordinates.resize(3, triDimensions[0]);
	    triMesh.Cell0DsId.reserve(triDimensions[0]);
	    triMesh.Cell1DsId.reserve(triDimensions[1]);
	    triMesh.Cell1DsExtrema.resize(2, triDimensions[1]);
	    triMesh.Cell2DsId.reserve(triDimensions[2]);
	    triMesh.Cell2DsEdges.reserve(triDimensions[2]);
	    triMesh.Cell2DsVertices.reserve(triDimensions[2]);
	
	    unsigned int vCount = 0, eCount = 0, fCount = 0;
	
	    std::map<std::pair<unsigned int, unsigned int>, unsigned int> edge_to_barycenter;
	    std::vector<unsigned int> barycenters;
	
	    for (unsigned int faceIdx = 0; faceIdx < baseMesh.Cell2DsId.size(); faceIdx++) {
	        const auto& faceVerts = baseMesh.Cell2DsVertices[faceIdx];
	        Vector3d A = baseMesh.Cell0DsCoordinates.col(faceVerts[0]);
	        Vector3d B = baseMesh.Cell0DsCoordinates.col(faceVerts[1]);
	        Vector3d C = baseMesh.Cell0DsCoordinates.col(faceVerts[2]);
	
	        std::vector<std::vector<unsigned int>> grid(level + 1);
	        for (unsigned int i = 0; i <= level; ++i) {
	            Vector3d from = ((double)i / level) * B + ((double)(level - i) / level) * A;
	            Vector3d to = ((double)i / level) * C + ((double)(level - i) / level) * A;
	            for (unsigned int j = 0; j <= i; ++j) {
	                Vector3d pos = (i == 0) ? A : ((double)j / i) * to + ((double)(i - j) / i) * from;
	                pos.normalize();
	                unsigned int id;
	                if (!VertexIsDupe(triMesh, pos, id)) {
	                    triMesh.Cell0DsId.push_back(vCount);
	                    triMesh.Cell0DsCoordinates.col(vCount) = pos;
	                    id = vCount++;
	                }
	                grid[i].push_back(id);
	            }
	        }
	
	        for (unsigned int i = 0; i < level; ++i) {
	            for (unsigned int j = 0; j <= i; ++j) {
	                // Triangle 1
	                std::array<unsigned int, 3> t1 = {grid[i][j], grid[i + 1][j], grid[i + 1][j + 1]};
	                Vector3d c = (triMesh.Cell0DsCoordinates.col(t1[0]) +
	                              triMesh.Cell0DsCoordinates.col(t1[1]) +
	                              triMesh.Cell0DsCoordinates.col(t1[2])) / 3.0;
	                c.normalize();
	                unsigned int c_id;
	                if (!VertexIsDupe(triMesh, c, c_id)) {
	                    triMesh.Cell0DsId.push_back(vCount);
	                    triMesh.Cell0DsCoordinates.col(vCount) = c;
	                    c_id = vCount++;
	                }
	                barycenters.push_back(c_id);
	                for (int k = 0; k < 3; ++k) {
	                    unsigned int from = t1[k];
	                    unsigned int to = t1[(k + 1) % 3];
	                    std::pair<unsigned int, unsigned int> edge = std::minmax(from, to);
	                    edge_to_barycenter[edge] = c_id;
	                    std::vector<unsigned int> verts = {from, to, c_id};
	                    triMesh.Cell2DsId.push_back(fCount++);
	                    triMesh.Cell2DsVertices.push_back(verts);
	                    std::vector<unsigned int> edges;
	                    for (int e = 0; e < 3; ++e) {
	                        unsigned int v1 = verts[e], v2 = verts[(e + 1) % 3];
	                        Vector2i extrema(v1, v2);
	                        if (v2 < v1) std::swap(extrema[0], extrema[1]);
	                        bool exists = false;
	                        for (unsigned int e_id = 0; e_id < eCount; ++e_id) {
	                            if ((triMesh.Cell1DsExtrema.col(e_id)(0) == extrema[0] &&
	                                 triMesh.Cell1DsExtrema.col(e_id)(1) == extrema[1])) {
	                                edges.push_back(e_id);
	                                exists = true;
	                                break;
	                            }
	                        }
	                        if (!exists) {
	                            triMesh.Cell1DsId.push_back(eCount);
	                            triMesh.Cell1DsExtrema.col(eCount) = extrema;
	                            edges.push_back(eCount++);
	                        }
	                    }
	                    triMesh.Cell2DsEdges.push_back(edges);
	                }
	                if (j < i) {
	                    // Triangle 2
	                    std::array<unsigned int, 3> t2 = {grid[i][j], grid[i][j + 1], grid[i + 1][j + 1]};
	                    c = (triMesh.Cell0DsCoordinates.col(t2[0]) +
	                         triMesh.Cell0DsCoordinates.col(t2[1]) +
	                         triMesh.Cell0DsCoordinates.col(t2[2])) / 3.0;
	                    c.normalize();
	                    if (!VertexIsDupe(triMesh, c, c_id)) {
	                        triMesh.Cell0DsId.push_back(vCount);
	                        triMesh.Cell0DsCoordinates.col(vCount) = c;
	                        c_id = vCount++;
	                    }
	                    barycenters.push_back(c_id);
	                    for (int k = 0; k < 3; ++k) {
	                        unsigned int from = t2[k];
	                        unsigned int to = t2[(k + 1) % 3];
	                        std::pair<unsigned int, unsigned int> edge = std::minmax(from, to);
	                        edge_to_barycenter[edge] = c_id;
	                        std::vector<unsigned int> verts = {from, to, c_id};
	                        triMesh.Cell2DsId.push_back(fCount++);
	                        triMesh.Cell2DsVertices.push_back(verts);
	                        std::vector<unsigned int> edges;
	                        for (int e = 0; e < 3; ++e) {
	                            unsigned int v1 = verts[e], v2 = verts[(e + 1) % 3];
	                            Vector2i extrema(v1, v2);
	                            if (v2 < v1) std::swap(extrema[0], extrema[1]);
	                            bool exists = false;
	                            for (unsigned int e_id = 0; e_id < eCount; ++e_id) {
	                                if ((triMesh.Cell1DsExtrema.col(e_id)(0) == extrema[0] &&
	                                     triMesh.Cell1DsExtrema.col(e_id)(1) == extrema[1])) {
	                                    edges.push_back(e_id);
	                                    exists = true;
	                                    break;
	                                }
	                            }
	                            if (!exists) {
	                                triMesh.Cell1DsId.push_back(eCount);
	                                triMesh.Cell1DsExtrema.col(eCount) = extrema;
	                                edges.push_back(eCount++);
	                            }
	                        }
	                        triMesh.Cell2DsEdges.push_back(edges);
	                    }
	                }
	            }
	        }
	    }
	
	    // Collega i baricentri adiacenti lungo ogni lato condiviso
	    for (const auto& [edge, b1] : edge_to_barycenter) {
	        auto range = edge_to_barycenter.equal_range(edge);
	        for (auto it = range.first; it != range.second; ++it) {
	            for (auto it2 = std::next(it); it2 != range.second; ++it2) {
	                unsigned int b2 = it2->second;
	                if (b1 != b2) {
	                    triMesh.Cell1DsId.push_back(eCount);
	                    triMesh.Cell1DsExtrema(0, eCount) = b1;
	                    triMesh.Cell1DsExtrema(1, eCount) = b2;
	                    triMesh.Cell2DsId.push_back(fCount);
	                    triMesh.Cell2DsVertices.push_back({b1, b2, edge.first});
	                    triMesh.Cell2DsEdges.push_back({eCount});
	                    eCount++;
	                    fCount++;
	                }
	            }
	        }
	    }
	
	    triMesh.Cell3DsId = {0};
	    triMesh.Cell3DsVertices = triMesh.Cell0DsId;
	    triMesh.Cell3DsEdges = triMesh.Cell1DsId;
	    triMesh.Cell3DsFaces = triMesh.Cell2DsId;
	    triMesh.NumCell0Ds = triMesh.Cell0DsId.size();
	    triMesh.NumCell1Ds = triMesh.Cell1DsId.size();
	    triMesh.NumCell2Ds = triMesh.Cell2DsId.size();
	    triMesh.NumCell3Ds = 1;
	
	    return true;
	}

}
