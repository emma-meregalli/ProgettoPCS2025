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

	// Funzione che controlla se un vertice è già nella lista
    bool VertexIsDupe(const PolyhedralMesh& mesh, const Vector3d& v, unsigned int& original_id){
        // Fisso una tolleranza per confrontare i vertici
	    double tol = 1e-12;
	    
        // Confronto con tutti i vertici già  inseriti nella lista
        for (size_t i = 0; i < mesh.Cell0DsId.size(); i++) {
            if ((mesh.Cell0DsCoordinates.col(i) - v).norm()/mesh.Cell0DsCoordinates.col(i).norm() < tol){    // Controlla l'errore relativo tra le norme dei vertici esistenti e quella del vertice nuovo
				original_id = i;    //Se il vertice esiste, allora restituisco il suo ID 
				return true;  
            }
        }
        return false;
    }

	// Funzione che controlla se uno spigolo è già nella lista
    bool EdgeIsDupe(const PolyhedralMesh& mesh, const Vector2i& e, unsigned int& original_id){
    	Vector2i v = e;
        for(size_t i = 0; i < mesh.Cell1DsId.size(); i++){
            if(mesh.Cell1DsExtrema.col(i) == v){
            	original_id = i;
                return true;
			}
            // Controllo i lati sia come (start, end) che come (end, start)
            swap(v[0],v[1]);
            if(mesh.Cell1DsExtrema.col(i) == v){
            	original_id = i;
                return true;
        	}
        }
        return false;
    }
    
    // Funzione di triangolazione - POLIEDRI GEODETICI DI CLASSE I
    bool GenerateTriangulatedMesh1(PolyhedralMesh& baseMesh, PolyhedralMesh& triMesh,      
        const unsigned int& b, const unsigned int& c, // Parametri della suddivisione
        const Vector3i& triDimensions) // Dimensione di (V,E,F) della mesh triangolata (con duplicati)
    {
        unsigned int level = b + c; // Numero di suddivisioni laterali per triangolo

        // Inizializzazione della struttura dati della mesh triangolata

        // Allocazione memoria per vertici (0D)
		triMesh.Cell0DsCoordinates = MatrixXd::Zero(3, triDimensions[0]); 
		triMesh.Cell0DsId.reserve(triDimensions[0]);
		
        // Allocazione memoria per lati (1D)
        triMesh.Cell1DsId.reserve(triDimensions[1]);
        triMesh.Cell1DsExtrema = MatrixXi::Zero(2, triDimensions[1]);
	
        // Allocazione memoria per facce (2D)
        triMesh.Cell2DsId.reserve(triDimensions[2]);
        triMesh.Cell2DsEdges.reserve(triDimensions[2]);
        triMesh.Cell2DsVertices.reserve(triDimensions[2]);
        for (auto& edgeList : triMesh.Cell2DsEdges){
            edgeList.resize(3); // Ogni faccia ha 3 spigoli
        }
		for (auto& vertList : triMesh.Cell2DsVertices){
        	vertList.resize(3); // Ogni faccia ha 3 vertici
		}	
		
        // Contatori per vertici, spigoli e facce
        unsigned int vCount = 0;
		unsigned int eCount = 0;
		unsigned int fCount = 0;

        // Triangolazione per ogni faccia del poliedro iniziale
        for (unsigned int faceIdx = 0; faceIdx < baseMesh.Cell2DsId.size(); faceIdx++) 
        {
            const auto& faceVerts = baseMesh.Cell2DsVertices[faceIdx]; //Prende i tre vertici della faccia corrente
        
            // Coordinate dei 3 vertici del triangolo originale
            Vector3d A = baseMesh.Cell0DsCoordinates.col(faceVerts[0]); // Vertice A
            Vector3d B = baseMesh.Cell0DsCoordinates.col(faceVerts[1]); // Vertice B
            Vector3d C = baseMesh.Cell0DsCoordinates.col(faceVerts[2]); // Vertice C
			
            vector<vector<unsigned int>> grid; // Griglia di vertici interni alla faccia
            // La griglia ha level = b + c righe e ogni riga i ha i + 1 elementi (forma triangolare)
            // Costruzione della griglia interpolata sulla faccia
            for (unsigned int i = 0; i <= level; i++){
                vector<unsigned int> row; // Riga corrente dei vertici 
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
						triMesh.Cell0DsId.push_back(vCount); // Memorizza l'ID del vertice
						for(unsigned int n = 0; n < 3; n++){
							triMesh.Cell0DsCoordinates(n,vCount) = pos(n); // Memorizza la posizione
						}
                    	row.push_back(vCount); // Aggiunge l'indice del vertice alla riga corrente
                    	vCount++; // Avanza il contatore dei vertici
					}
					else{
						row.push_back(original_id); // Se il vertice esiste già
					}
                }
                grid.push_back(row); // Aggiunge la riga corrente alla griglia
            }

            //Crea i nuovi lati dati dalla triangolazione e aggiorna la lista dei lati
            unsigned int original_id = 0;
            
            for(size_t i = 0; i < grid.size(); i++){
                Vector2i extrema;
                for(size_t j = 0; j <= i; j++){
                    if(i < grid.size() - 1){
                    	// Lato sotto a sinistra
                    	extrema[0] = grid[i][j];
						extrema[1] = grid[i + 1][j];
                    	
                    	// Aggiorna la lista dei lati, se il vertice non esiste già
                    	if(!EdgeIsDupe(triMesh, extrema, original_id)){
                    		triMesh.Cell1DsId.push_back(eCount);
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i + 1][j]; 
	                        eCount++;	
						}
						
						// Lato sotto a destra
						extrema[0] = grid[i][j];
						extrema[1] = grid[i + 1][j+1];
						
						// Aggiorna la lista dei lati, se il vertice non esiste già
						if(!EdgeIsDupe(triMesh, extrema, original_id)){
	                        triMesh.Cell1DsId.push_back(eCount);
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i + 1][j + 1];
	                        eCount++;
	                    }
                    }
                    
                    // Lato accanto a destra
                    if(i > 0 && j <= i - 1){
                    	extrema[0] = grid[i][j];
						extrema[1] = grid[i][j + 1];
						
						// Aggiorna la lista dei lati, se il vertice non esiste già
                    	if(!EdgeIsDupe(triMesh, extrema, original_id)){
	                        triMesh.Cell1DsId.push_back(eCount);
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i][j + 1];
	                        eCount++;
                    	}
                    }
                }
            }

			//Crea le nuove facce dopo la triangolazione
			for (size_t i = 0; i < grid.size() - 1; i++){
				for (size_t j = 0; j <= i; j++){
					
					// Facce con la punta verso l'alto
					vector<unsigned int> v1 = {grid[i][j],grid[i + 1][j],grid[i + 1][j + 1]}; 
					vector<unsigned int> e1;
					// Memorizza i lati corrispodenti
					for(unsigned int v = 0; v < 3; v++){
						unsigned int from = v1[v];
						unsigned int to = v1[(v + 1) % 3];
						for(unsigned int k = 0; k < triMesh.Cell1DsId.size(); k++){
							if((from == triMesh.Cell1DsExtrema(0, k) && to == triMesh.Cell1DsExtrema(1, k)) || (from == triMesh.Cell1DsExtrema(1, k) && to == triMesh.Cell1DsExtrema(0, k))){
								e1.push_back(k);
								break;
							}
						}
					}
					// Aggiorna le variabili relative a Cell2Ds
					triMesh.Cell2DsId.push_back(fCount);
					triMesh.Cell2DsVertices.push_back(v1);
					triMesh.Cell2DsEdges.push_back(e1);
					fCount++;
					
					// Facce con la punta verso il basso
					if(i > 0 && j <= i - 1){
						vector<unsigned int> v2 = {grid[i][j],grid[i][j + 1],grid[i + 1][j + 1]};
						vector<unsigned int> e2;
						// Memorizza i lati corrispodenti
						for(unsigned int v = 0; v < 3; v++){
							unsigned int from = v2[v];
							unsigned int to = v2[(v + 1) % 3];
							for(unsigned int k = 0; k < triMesh.Cell1DsId.size(); k++){
								if((from == triMesh.Cell1DsExtrema(0, k) && to == triMesh.Cell1DsExtrema(1, k)) || (from == triMesh.Cell1DsExtrema(1, k) && to == triMesh.Cell1DsExtrema(0, k))){
									e2.push_back(k);
									break;
								}
							}	
						}
						// Aggiorna le variabili relative a Cell2Ds
						triMesh.Cell2DsId.push_back(fCount);
						triMesh.Cell2DsVertices.push_back(v2);
						triMesh.Cell2DsEdges.push_back(e2);
						fCount++;
					}		
			    }
		    }
		}
    
		//Crea il poliedro triangolato
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
	
	// Funzione di triangolazione - POLIEDRI GEODETICI DI CLASSE II
    bool GenerateTriangulatedMesh2(PolyhedralMesh& baseMesh, PolyhedralMesh& triMesh,
        const unsigned int& b, const unsigned int& c, // Parametri della suddivisione
        const Vector3i& triDimensions) // Dimensione di (V,E,F) della mesh triangolata
    {
        unsigned int level = b; // Numero di suddivisioni laterali per triangolo (per la griglia base)

        // Inizializzazione della struttura dati della mesh triangolata
        triMesh.Cell0DsCoordinates.resize(3, triDimensions[0]);
        triMesh.Cell0DsId.reserve(triDimensions[0]);
        
        triMesh.Cell1DsId.reserve(triDimensions[1]);
        triMesh.Cell1DsExtrema.resize(2, triDimensions[1]);
        
        triMesh.Cell2DsId.reserve(triDimensions[2]);
        triMesh.Cell2DsEdges.reserve(triDimensions[2]);
        triMesh.Cell2DsVertices.reserve(triDimensions[2]);
        
        // Contatori per vertici, spigoli e facce
        unsigned int vCount = 0;
        unsigned int eCount = 0;
        unsigned int fCount = 0;
        
        // Triangolazione per ogni faccia del poliedro iniziale
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
                for (unsigned int j = 0; j <= i; j++){ 
                    Vector3d pos;
                    if (i == 0){
                        pos = A;
                    } 
					else{
                        pos = ((double)j / i) * to + ((double)(i - j) / i) * from;
                    }
                    //pos = pos / pos.norm();
                    if(!VertexIsDupe(triMesh, pos, original_id)){
                        triMesh.Cell0DsId.push_back(vCount);
                        triMesh.Cell0DsCoordinates.col(vCount) = pos;
                        row.push_back(vCount);
                        vCount++;
                    } 
					else{
                        row.push_back(original_id);
                    }
                }
                grid_base_verts.push_back(row);
            }
            
            vector<unsigned int> barycenters;  // Vettore per segnare i baricentri dei triangoli verso l'alto della riga corrente
            vector<unsigned int> barycenters2;  // Vettore per segnare i baricentri dei triangoli verso il basso della riga corrente
        	vector<vector<unsigned int>> barycenters_grid; // Griglia dei baricentri
        	vector<vector<unsigned int>> barycenters_grid2; // Griglia dei baricentri 2

            // Ora per ogni triangolo creato con la triangolazione 1 applico la triangolazione 2
            for (unsigned int i = 0; i <= level; i++) {
                for (unsigned int j = 0; j <= i; j++) {
                	
                    // Triangoli con la punto verso l'alto
                    if(i != level) {
                    	vector<unsigned int> triangleVertices1 = {grid_base_verts[i][j], grid_base_verts[i + 1][j], grid_base_verts[i + 1][j + 1]};
                    
	                    // Prende le coordinate della faccia corrente ottenuta dalla triangolazione 1
	                    Vector3d p1_coord = triMesh.Cell0DsCoordinates.col(triangleVertices1[0]);
	                    Vector3d p2_coord = triMesh.Cell0DsCoordinates.col(triangleVertices1[1]);
	                    Vector3d p3_coord = triMesh.Cell0DsCoordinates.col(triangleVertices1[2]);
	                    
	                    // Per ogni caso bisogna vedere se i lati 12, 23, 31 vanno tenuti oppure no
	                    bool exists12 = false;
	                    bool exists23 = false;
	                    bool exists31 = false;
	
	                    // Calcola i punti medi di ogni lato (se vanno tenuti) e il baricentro del triangolo
	                    Vector3d mid12_pos, mid31_pos, mid23_pos;
	                    
	                    if(j == 0) {
	                    	mid12_pos = (p1_coord + p2_coord) / 2.0;
	                    	//mid12_pos =  mid12_pos / mid12_pos.norm();
	                    	exists12 = true;
						}
						if(j == i) {
							mid31_pos = (p3_coord + p1_coord) / 2.0;
							//mid31_pos =  mid31_pos / mid31_pos.norm();
							exists31 = true;
						}
						if(i == level - 1){
							mid23_pos = (p2_coord + p3_coord) / 2.0;
							//mid23_pos =  mid23_pos / mid23_pos.norm();
							exists23 = true;
						}
	                    
	                    Vector3d barycenter_pos = (p1_coord + p2_coord + p3_coord) / 3.0;
						//barycenter_pos = barycenter_pos / barycenter_pos.norm();
						
	                    unsigned int mid12_id, mid23_id, mid31_id, barycenter_id;
	                    unsigned int original_id;
						
	                    // Aggiunge i punti medi se non esistono già
	                    if(exists12){
	                    	if(!VertexIsDupe(triMesh, mid12_pos, original_id)){
		                        triMesh.Cell0DsId.push_back(vCount);
		                        triMesh.Cell0DsCoordinates.col(vCount) = mid12_pos;
		                        mid12_id = vCount;
		                        vCount++;
	                    	} 
							else{ 
								mid12_id = original_id; 
							}
						}
						
						if(exists23){
		                    if(!VertexIsDupe(triMesh, mid23_pos, original_id)){
		                        triMesh.Cell0DsId.push_back(vCount);
		                        triMesh.Cell0DsCoordinates.col(vCount) = mid23_pos;
		                        mid23_id = vCount;
		                        vCount++;
		                    } 
							else{
								mid23_id = original_id; 
							}
						}
						
						if(exists31){
		                    if(!VertexIsDupe(triMesh, mid31_pos, original_id)){
		                        triMesh.Cell0DsId.push_back(vCount);
		                        triMesh.Cell0DsCoordinates.col(vCount) = mid31_pos;
		                        mid31_id = vCount;
		                        vCount++;
		                    } 
							else{ 
								mid31_id = original_id; 
							}
						}
	                    
	                    // I baricentri non rischiano di essere copie
	                    triMesh.Cell0DsId.push_back(vCount);
	                    triMesh.Cell0DsCoordinates.col(vCount) = barycenter_pos;
	                    barycenter_id = vCount;
	                    barycenters.push_back(vCount);
	                    vCount++;
	                    
	                    // Crea i nuovi triangoli vertice-puntomedio-baricentro (se sono previsti per questo triangolo)
	                    vector<vector<unsigned int>> new_sub_triangles_1;
						
						if(exists12){
							new_sub_triangles_1.push_back(std::vector<unsigned int>{triangleVertices1[0], mid12_id, barycenter_id});
							new_sub_triangles_1.push_back(std::vector<unsigned int>{mid12_id, triangleVertices1[1], barycenter_id});
						}
						if(exists23){
							new_sub_triangles_1.push_back(std::vector<unsigned int>{triangleVertices1[1], mid23_id, barycenter_id});
							new_sub_triangles_1.push_back(std::vector<unsigned int>{mid23_id, triangleVertices1[2], barycenter_id});
						}
						if(exists31){
							new_sub_triangles_1.push_back(std::vector<unsigned int>{triangleVertices1[2], mid31_id, barycenter_id});
							new_sub_triangles_1.push_back(std::vector<unsigned int>{mid31_id, triangleVertices1[0], barycenter_id});
						}
						
						
						// Se il triangolo è l'ultimo del rispettivo strato, aggiungo i triangoli che si creano collegando il suo baricentro a quello adiacente a sinistra
						if(i > 0 && j == i){
							new_sub_triangles_1.push_back(std::vector<unsigned int>{barycenters2.back(), barycenter_id, triangleVertices1[0]});
							new_sub_triangles_1.push_back(std::vector<unsigned int>{barycenters2.back(), barycenter_id, triangleVertices1[1]});
						}
	
						unsigned int original_id2;
						
						// Controlla se i lati esistono già prima di aggiungerli
	                    for(const auto& new_verts : new_sub_triangles_1){
	                        vector<unsigned int> current_edges;
	                        for(unsigned int k = 0; k < 3; k++){
	                            Vector2i edge_extrema;
	                            edge_extrema[0] = new_verts[k];
	                            edge_extrema[1] = new_verts[(k + 1) % 3];
	                            
	                            if(EdgeIsDupe(triMesh, edge_extrema, original_id2)){
	                            	current_edges.push_back(original_id2);
								}
	                            else {
	                                triMesh.Cell1DsId.push_back(eCount);
	                                triMesh.Cell1DsExtrema.col(eCount) = edge_extrema;
	                                current_edges.push_back(eCount);
	                                eCount++;
	                            }
	                        }
	                        
	                        // Aggiorna Cell2Ds
	                        triMesh.Cell2DsId.push_back(fCount);
	                        triMesh.Cell2DsVertices.push_back(new_verts);
	                        triMesh.Cell2DsEdges.push_back(current_edges);
	                        fCount++;
	                    }
	            	}
                    
                    // Triangoli con la punta verso il basso (si parte dal secondo strato)
                    if (i > 0 && j < i) {
                        vector<unsigned int> triangleVertices2 = {grid_base_verts[i][j], grid_base_verts[i][j + 1], grid_base_verts[i + 1][j + 1]};
                        
                        // Prende le coordinate della faccia corrente ottenuta dalla triangolazione 1
                        Vector3d p1_coord = triMesh.Cell0DsCoordinates.col(triangleVertices2[0]);
                        Vector3d p2_coord = triMesh.Cell0DsCoordinates.col(triangleVertices2[1]);
                        Vector3d p3_coord = triMesh.Cell0DsCoordinates.col(triangleVertices2[2]);
						
						// In questo caso salva solo le coordinate del baricentro
                        Vector3d barycenter_pos = (p1_coord + p2_coord + p3_coord) / 3.0;
                        //barycenter_pos = barycenter_pos / barycenter_pos.norm();

                        triMesh.Cell0DsId.push_back(vCount);
	                    triMesh.Cell0DsCoordinates.col(vCount) = barycenter_pos;
	                    unsigned int barycenter_id = vCount;
	                    barycenters2.push_back(vCount);
	                    vCount++;

						// Aggiunge i triangoli che si creano dal collegamento col baricentro del triangolo sopra di quello a sinistra
                        vector<vector<unsigned int>> new_sub_triangles_2 = {
                            {triangleVertices2[0], barycenters_grid[i - 1][j], barycenter_id},
                        	{triangleVertices2[1], barycenters_grid[i - 1][j], barycenter_id},
                        	{triangleVertices2[0], barycenters.back(), barycenter_id},
                        	{triangleVertices2[2], barycenters.back(), barycenter_id}
                        };
                       
                       unsigned int original_id2;
                       
                       // Controlla se i lati esistono già prima di aggiungerli
                        for(const auto& new_verts : new_sub_triangles_2){
                            vector<unsigned int> current_edges;
                            for(unsigned int k = 0; k < 3; k++){
                                Vector2i edge_extrema;
                                edge_extrema[0] = new_verts[k];
                                edge_extrema[1] = new_verts[(k + 1) % 3];
                                
                                if(EdgeIsDupe(triMesh, edge_extrema, original_id2)){
	                            	current_edges.push_back(original_id2);
								}
	                            else {
	                                triMesh.Cell1DsId.push_back(eCount);
	                                triMesh.Cell1DsExtrema.col(eCount) = edge_extrema;
	                                current_edges.push_back(eCount);
	                                eCount++;
	                            }
                            }
                            
                            // Aggiorna Cell2Ds
                            triMesh.Cell2DsId.push_back(fCount);
	                        triMesh.Cell2DsVertices.push_back(new_verts);
	                        triMesh.Cell2DsEdges.push_back(current_edges);
	                        fCount++;
                        }  
                    }
                }
                // Aggiorna la griglia dei baricentri, necessaria per creare gli spigoli tra di essi
                barycenters_grid.push_back(barycenters);
                barycenters_grid2.push_back(barycenters2);
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
    	}     	
		return true;
	}	
}
