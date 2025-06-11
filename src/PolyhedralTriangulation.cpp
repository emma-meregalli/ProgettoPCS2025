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
                    // pos = pos/pos.norm();
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
        unsigned int level = b + c; // Numero di suddivisioni laterali per triangolo (per la griglia base)

        // Inizializzazione della struttura dati della mesh triangolata
        triMesh.Cell0DsCoordinates.resize(3, triDimensions[0]);
        triMesh.Cell0DsId.reserve(triDimensions[0]);
        triMesh.Cell1DsId.reserve(triDimensions[1]);
        triMesh.Cell1DsExtrema.resize(2, triDimensions[1]);
        triMesh.Cell2DsId.reserve(triDimensions[2]);
        triMesh.Cell2DsEdges.reserve(triDimensions[2]);
        triMesh.Cell2DsVertices.reserve(triDimensions[2]);
        
        unsigned int vCount = 0;
        unsigned int eCount = 0;
        unsigned int fCount = 0;
        
        for (unsigned int faceIdx = 0; faceIdx < baseMesh.Cell2DsId.size(); faceIdx++) 
        {
            const auto& faceVerts = baseMesh.Cell2DsVertices[faceIdx];
            Vector3d A = baseMesh.Cell0DsCoordinates.col(faceVerts[0]);
            Vector3d B = baseMesh.Cell0DsCoordinates.col(faceVerts[1]);
            Vector3d C = baseMesh.Cell0DsCoordinates.col(faceVerts[2]);
            
            // Griglia di vertici base (come in GenerateTriangulatedMesh1)
            vector<vector<unsigned int>> grid_base_verts; 
            for (unsigned int i = 0; i <= level; i++){
                vector<unsigned int> row;
                unsigned int original_id;
                Vector3d from = ((double)i / level) * B + ((double)(level - i) / level) * A;
                Vector3d to = ((double)i / level) * C + ((double)(level - i) / level) * A;
                for (unsigned int j = 0; j <= i; j++) { 
                    Vector3d pos;
                    if (i == 0) {
                        pos = A;
                    } else {
                        pos = ((double)j / i) * to + ((double)(i - j) / i) * from;
                    }
                    
                    if(!VertexIsDupe(triMesh, pos, original_id)){
                        triMesh.Cell0DsId.push_back(vCount);
                        triMesh.Cell0DsCoordinates.col(vCount) = pos;
                        row.push_back(vCount);
                        vCount++;
                    } else {
                        row.push_back(original_id);
                    }
                }
                grid_base_verts.push_back(row);
            }

            // Ora, iteriamo sui triangoli creati dalla griglia base (come se fosse la GenerateTriangulatedMesh1)
            // e per ognuno di essi, applichiamo la raffinazione di Classe II.
            for (unsigned int i = 0; i < level; ++i) {
                for (unsigned int j = 0; j <= i; ++j) {
					cout<<vCount<<", "<<eCount<<", "<<fCount<<endl;
                    // Raffinazione per il Triangolo "verso il basso"
                    vector<unsigned int> current_triangle_verts_1 = {grid_base_verts[i][j], grid_base_verts[i+1][j], grid_base_verts[i+1][j+1]};
                    
                    // Get coordinates of the base triangle vertices
                    Vector3d p1_coord = triMesh.Cell0DsCoordinates.col(current_triangle_verts_1[0]);
                    Vector3d p2_coord = triMesh.Cell0DsCoordinates.col(current_triangle_verts_1[1]);
                    Vector3d p3_coord = triMesh.Cell0DsCoordinates.col(current_triangle_verts_1[2]);

                    // Calculate midpoints and barycenter for this triangle
                    Vector3d mid12_pos = (p1_coord + p2_coord) / 2.0;
                    Vector3d mid23_pos = (p2_coord + p3_coord) / 2.0;
                    Vector3d mid31_pos = (p3_coord + p1_coord) / 2.0;
                    Vector3d barycenter_pos = (p1_coord + p2_coord + p3_coord) / 3.0;

                    unsigned int mid12_id, mid23_id, mid31_id, barycenter_id;
                    unsigned int original_id;

                    // Add midpoints and barycenter to triMesh if not duplicates
                    if(!VertexIsDupe(triMesh, mid12_pos, original_id)){
                        triMesh.Cell0DsId.push_back(vCount);
                        triMesh.Cell0DsCoordinates.col(vCount) = mid12_pos;
                        mid12_id = vCount;
                        vCount++;
                    } else { mid12_id = original_id; }

                    if(!VertexIsDupe(triMesh, mid23_pos, original_id)){
                        triMesh.Cell0DsId.push_back(vCount);
                        triMesh.Cell0DsCoordinates.col(vCount) = mid23_pos;
                        mid23_id = vCount;
                        vCount++;
                    } else { mid23_id = original_id; }

                    if(!VertexIsDupe(triMesh, mid31_pos, original_id)){
                        triMesh.Cell0DsId.push_back(vCount);
                        triMesh.Cell0DsCoordinates.col(vCount) = mid31_pos;
                        mid31_id = vCount;
                        vCount++;
                    } else { mid31_id = original_id; }
                    
                    if(!VertexIsDupe(triMesh, barycenter_pos, original_id)){
                        triMesh.Cell0DsId.push_back(vCount);
                        triMesh.Cell0DsCoordinates.col(vCount) = barycenter_pos;
                        barycenter_id = vCount;
                        vCount++;
                    } else { barycenter_id = original_id; }
                    
                    // Create the new 6 triangles for current_triangle_verts_1
                    vector<vector<unsigned int>> new_sub_triangles_1 = {
                        {current_triangle_verts_1[0], mid12_id, barycenter_id},
                        {mid12_id, current_triangle_verts_1[1], barycenter_id},
                        {current_triangle_verts_1[1], mid23_id, barycenter_id},
                        {mid23_id, current_triangle_verts_1[2], barycenter_id},
                        {current_triangle_verts_1[2], mid31_id, barycenter_id},
                        {mid31_id, current_triangle_verts_1[0], barycenter_id}
                    };

                    for(const auto& new_verts : new_sub_triangles_1){
                        vector<unsigned int> current_edges;
                        for(unsigned int k=0; k<3; ++k){
                            Vector2i edge_extrema;
                            edge_extrema[0] = new_verts[k];
                            edge_extrema[1] = new_verts[(k+1)%3];
                            
                            bool found_edge = false;
                            for(unsigned int existing_e_id = 0; existing_e_id < triMesh.Cell1DsId.size(); ++existing_e_id){
                                if((triMesh.Cell1DsExtrema(0, existing_e_id) == edge_extrema[0] && triMesh.Cell1DsExtrema(1, existing_e_id) == edge_extrema[1]) ||
                                   (triMesh.Cell1DsExtrema(0, existing_e_id) == edge_extrema[1] && triMesh.Cell1DsExtrema(1, existing_e_id) == edge_extrema[0])){
                                    current_edges.push_back(triMesh.Cell1DsId[existing_e_id]);
                                    found_edge = true;
                                    break;
                                }
                            }
                            if (!found_edge) {
                                triMesh.Cell1DsId.push_back(eCount);
                                triMesh.Cell1DsExtrema.col(eCount) = edge_extrema;
                                current_edges.push_back(eCount);
                                eCount++;
                            }
                        }
                        triMesh.Cell2DsId.push_back(fCount);
                        triMesh.Cell2DsVertices.push_back(new_verts);
                        triMesh.Cell2DsEdges.push_back(current_edges);
                        fCount++;
                    }

                    // Handle Triangolo "verso l'alto" (if it exists)
                    if (j < i) {
                        vector<unsigned int> current_triangle_verts_2 = {grid_base_verts[i][j], grid_base_verts[i][j+1], grid_base_verts[i+1][j+1]};
                        
                        p1_coord = triMesh.Cell0DsCoordinates.col(current_triangle_verts_2[0]);
                        p2_coord = triMesh.Cell0DsCoordinates.col(current_triangle_verts_2[1]);
                        p3_coord = triMesh.Cell0DsCoordinates.col(current_triangle_verts_2[2]);

                        mid12_pos = (p1_coord + p2_coord) / 2.0;
                        mid23_pos = (p2_coord + p3_coord) / 2.0;
                        mid31_pos = (p3_coord + p1_coord) / 2.0;
                        barycenter_pos = (p1_coord + p2_coord + p3_coord) / 3.0;

                        // Add midpoints and barycenter to triMesh if not duplicates
                        if(!VertexIsDupe(triMesh, mid12_pos, original_id)){
                            triMesh.Cell0DsId.push_back(vCount);
                            triMesh.Cell0DsCoordinates.col(vCount) = mid12_pos;
                            mid12_id = vCount;
                            vCount++;
                        } else { mid12_id = original_id; }

                        if(!VertexIsDupe(triMesh, mid23_pos, original_id)){
                            triMesh.Cell0DsId.push_back(vCount);
                            triMesh.Cell0DsCoordinates.col(vCount) = mid23_pos;
                            mid23_id = vCount;
                            vCount++;
                        } else { mid23_id = original_id; }

                        if(!VertexIsDupe(triMesh, mid31_pos, original_id)){
                            triMesh.Cell0DsId.push_back(vCount);
                            triMesh.Cell0DsCoordinates.col(vCount) = mid31_pos;
                            mid31_id = vCount;
                            vCount++;
                        } else { mid31_id = original_id; }
                        
                        if(!VertexIsDupe(triMesh, barycenter_pos, original_id)){
                            triMesh.Cell0DsId.push_back(vCount);
                            triMesh.Cell0DsCoordinates.col(vCount) = barycenter_pos;
                            barycenter_id = vCount;
                            vCount++;
                        } else { barycenter_id = original_id; }

                        vector<vector<unsigned int>> new_sub_triangles_2 = {
                            {current_triangle_verts_2[0], mid12_id, barycenter_id},
                            {mid12_id, current_triangle_verts_2[1], barycenter_id},
                            {current_triangle_verts_2[1], mid23_id, barycenter_id},
                            {mid23_id, current_triangle_verts_2[2], barycenter_id},
                            {current_triangle_verts_2[2], mid31_id, barycenter_id},
                            {mid31_id, current_triangle_verts_2[0], barycenter_id}
                        };
                        
                        for(const auto& new_verts : new_sub_triangles_2){
                            vector<unsigned int> current_edges;
                            for(unsigned int k=0; k<3; ++k){
                                Vector2i edge_extrema;
                                edge_extrema[0] = new_verts[k];
                                edge_extrema[1] = new_verts[(k+1)%3];
                                
                                bool found_edge = false;
                                for(unsigned int existing_e_id = 0; existing_e_id < triMesh.Cell1DsId.size(); ++existing_e_id){
                                    if((triMesh.Cell1DsExtrema(0, existing_e_id) == edge_extrema[0] && triMesh.Cell1DsExtrema(1, existing_e_id) == edge_extrema[1]) ||
                                       (triMesh.Cell1DsExtrema(0, existing_e_id) == edge_extrema[1] && triMesh.Cell1DsExtrema(1, existing_e_id) == edge_extrema[0])){
                                        current_edges.push_back(triMesh.Cell1DsId[existing_e_id]);
                                        found_edge = true;
                                        break;
                                    }
                                }
                                if (!found_edge) {
                                    triMesh.Cell1DsId.push_back(eCount);
                                    triMesh.Cell1DsExtrema.col(eCount) = edge_extrema;
                                    current_edges.push_back(eCount);
                                    eCount++;
                                }
                            }
                            triMesh.Cell2DsId.push_back(fCount);
                            triMesh.Cell2DsVertices.push_back(new_verts);
                            triMesh.Cell2DsEdges.push_back(current_edges);
                            fCount++;
                        }
                    }
                }
            }
        }
    
        // Aggiorna Cell3DsId, Cell3DsVertices, Cell3DsEdges, Cell3DsFaces
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
