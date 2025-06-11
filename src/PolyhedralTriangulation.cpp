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
        // Confronto con tutti i vertici gi√† inseriti nella lista
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
                vector<unsigned int> bars; //riga corrente dei baricentri
				unsigned int original_id;
                // Calcolo il punto iniziale e finale della riga i-esima 
                // Partiziono il lato in base al valore di b e c
                Vector3d from = ((double)i / level) * B + ((double)(level - i) / level) * A;
                Vector3d to = ((double)i / level) * C + ((double)(level - i) / level) * A;
                for (unsigned int j = 0; j <= i; j++) { 
                    // Interpolo tra from e to per ottenere un punto interno
                    Vector3d pos;
                    Vector3d bar;
                    if (i == 0) {
                        pos = A;
                    } else {
                        pos = ((double)j / i) * to + ((double)(i - j) / i) * from;
                    }
                    // pos = pos/pos.norm();
                    if(i!=0){
                    	if(j==0){
                    		bar = (row[i-1][0]+pos)/2;
						}
						if(j==i){
                    		bar = (row[i-1][i-1]+pos)/2;
						}
						
						//Controllo duplicati dei baricentri dei nuovi lati
						if(!VertexIsDupe(triMesh, bar, original_id)){
						triMesh.Cell0DsId.push_back(vCount);           // Salva ID
							for(unsigned int n=0; n<3; n++){
								triMesh.Cell0DsCoordinates(n,vCount) = bar(n); // Salva posizione
							}
	                    	bars.push_back(vCount); // Aggiungi indice del vertice alla riga corrente
	                    	vCount++; // Avanza contatore vertice
						}
						else{
							bars.push_back(original_id);
						}
					}
					
					if(i==level){
						
					}
					
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
                grid.push_back(bars); // Aggiungi riga alla griglia
                grid.push_back(row); // Aggiungi riga alla griglia
			}
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
        unsigned int level = b; // Numero di livelli della faccia triangolata

        // Inizializzazione della struttura dati della mesh triangolata (identica a GenerateTriangulatedMesh1)
        triMesh.Cell0DsCoordinates = MatrixXd::Zero(3,triDimensions[0]);
        triMesh.Cell0DsId.reserve(triDimensions[0]);
        
        triMesh.Cell1DsId.reserve(triDimensions[1]);
        triMesh.Cell1DsExtrema = MatrixXi::Zero(2, triDimensions[1]);
        
        triMesh.Cell2DsId.reserve(triDimensions[2]);
        triMesh.Cell2DsEdges.reserve(triDimensions[2]);
        triMesh.Cell2DsVertices.reserve(triDimensions[2]);
        
        unsigned int vCount = 0;
        unsigned int eCount = 0;
        unsigned int fCount = 0;
        
        // Ciclo su tutte le facce della mesh di base
        for (unsigned int faceIdx = 0; faceIdx < baseMesh.Cell2DsId.size(); faceIdx++)
        {
            const auto& faceVerts = baseMesh.Cell2DsVertices[faceIdx];
            
            // Coordinate dei 3 vertici del triangolo originale
            Vector3d A = baseMesh.Cell0DsCoordinates.col(faceVerts[0]); // Vertice A
            Vector3d B = baseMesh.Cell0DsCoordinates.col(faceVerts[1]); // Vertice B
            Vector3d C = baseMesh.Cell0DsCoordinates.col(faceVerts[2]); // Vertice C
            
            vector<vector<unsigned int>> grid; // Griglia di vertici interni alla faccia
            
            // Costruzione della griglia interpolata sulla faccia per Classe II
            // P_ij = (b-i-j)/b * A + i/b * B + j/b * C
            for (unsigned int i = 0; i <= level; i++){ // i goes from 0 to b
                vector<unsigned int> row;
                unsigned int original_id;
                for (unsigned int j = 0; j <= level - i; ++j) { // j goes from 0 to b-i
                    unsigned int k = level - i - j; // k is the third barycentric coordinate

                    Vector3d pos = ((double)k / level) * A + ((double)i / level) * B + ((double)j / level) * C;
                    pos = pos / pos.norm(); // Normalizza sul raggio unitario

                    if(!VertexIsDupe(triMesh, pos, original_id)){
                        triMesh.Cell0DsId.push_back(vCount);
                        triMesh.Cell0DsCoordinates.col(vCount) = pos; // Assign entire column
                        row.push_back(vCount);
                        vCount++;
                    }
                    else{
                        row.push_back(original_id);
                    }
                }
                grid.push_back(row);
            }

            // Creiamo i nuovi lati e le nuove facce basandoci sulla griglia
            // La logica di creazione di lati e facce Ë leggermente pi˘ complessa per la Classe II
            // perchÈ la griglia interna non Ë un semplice "quadrato" come nella Classe I.
            // Invece di due triangoli per cella, qui si considerano triangoli orientati in modo diverso.
            // I loop per creare triangoli saranno simili, ma gli indici della 'grid' cambieranno.

            // First, add all edges for the current face's grid
            for (size_t i = 0; i < grid.size(); ++i) { // i is the row index
                for (size_t j = 0; j < grid[i].size(); ++j) { // j is the column index within the row
                    Vector2i extrema;
                    unsigned int v_curr = grid[i][j];

                    // Edge to the right (if not last element in row)
                    if (j + 1 < grid[i].size()) {
                        extrema[0] = v_curr;
                        extrema[1] = grid[i][j+1];
                        if (!EdgeIsDupe(triMesh, extrema)) {
                            triMesh.Cell1DsId.push_back(eCount);
                            triMesh.Cell1DsExtrema.col(eCount) = extrema;
                            eCount++;
                        }
                    }

                    // Edge diagonally down-left (if not last row)
                    if (i + 1 < grid.size()) {
                        extrema[0] = v_curr;
                        extrema[1] = grid[i+1][j]; // This vertex is in the next row, same 'j' coordinate for Class II grid
                        if (!EdgeIsDupe(triMesh, extrema)) {
                            triMesh.Cell1DsId.push_back(eCount);
                            triMesh.Cell1DsExtrema.col(eCount) = extrema;
                            eCount++;
                        }
                    }
                }
            }


            // Now, create the new faces for the current base face
            for (size_t i = 0; i < grid.size() - 1; ++i) {
                for (size_t j = 0; j < grid[i].size(); ++j) { // Note: j < grid[i].size()
                    // Triangle 1: pointing "up" (relative to the bottom edge of the current "quad")
                    // Vertices: grid[i][j], grid[i+1][j], grid[i+1][j+1]
                    // (This forms the "bottom-left" triangle of a square-like region if j is not at end)
                    if (j + 1 < grid[i+1].size()) { // Ensure grid[i+1][j+1] exists
                        vector<unsigned int> v1 = {grid[i][j], grid[i+1][j], grid[i+1][j+1]};
                        vector<unsigned int> e1;
                        for(unsigned int v_idx=0; v_idx<3; v_idx++){
                            unsigned int from=v1[v_idx];
                            unsigned int to= v1[(v_idx + 1) % 3];
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
                    }

                    // Triangle 2: pointing "down"
                    // Vertices: grid[i][j], grid[i][j+1], grid[i+1][j+1]
                    // (This forms the "top-right" triangle of a square-like region if j is not at end)
                    // This is for the "diamond" pattern, a bit different from Class I
                    if (j + 1 < grid[i].size() && j + 1 < grid[i+1].size()) { // Ensure grid[i][j+1] and grid[i+1][j+1] exist
                        vector<unsigned int> v2 = {grid[i][j], grid[i][j+1], grid[i+1][j+1]};
                        vector<unsigned int> e2;
                        for(unsigned int v_idx=0; v_idx<3; v_idx++){
                            unsigned int from=v2[v_idx];
                            unsigned int to=v2[(v_idx+1)%3];
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
        
        triMesh.NumCell0Ds = triMesh.Cell0DsId.size();
        triMesh.NumCell1Ds = triMesh.Cell1DsId.size();
        triMesh.NumCell2Ds = triMesh.Cell2DsId.size();
        triMesh.NumCell3Ds = 1;

        return true;
    }

   
}
