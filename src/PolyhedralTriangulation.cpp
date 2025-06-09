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

    bool VertexIsDupe(const PolyhedralMesh& mesh, const Vector3d& v){
        //Fisso una tolleranza per confrontare i vertici
	    double tol=1e-12;
        // Confronto con tutti i vertici già inseriti nella lista
        for (size_t i = 0; i < mesh.Cell0DsId.size(); i++) {
            if ((mesh.Cell0DsCoordinates.col(i) - v).norm() < tol){    //Se il vertice esiste, allora restituisco il suo ID (bisogna fare un controllo con la tolleranza?)
				return true;  
            }
        }
        return false;
    }

    bool EdgeIsDupe(const PolyhedralMesh& mesh, Vector2i& e){
        for(size_t i=0; i<mesh.Cell1DsId.size(); i++){
            if(mesh.Cell1DsExtrema.col(i)==e)
                return true;
            swap(e[0],e[1]);
            if(mesh.Cell1DsExtrema.col(i)==e)
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
        triMesh.Cell2DsId.resize(triDimensions[2]);
        triMesh.Cell2DsEdges.resize(triDimensions[2]);
        triMesh.Cell2DsVertices.resize(triDimensions[2]);
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
					if(!VertexIsDupe(triMesh, pos)){
						triMesh.Cell0DsCoordinates.col(vCount) = pos; // Salva posizione
                    	triMesh.Cell0DsId.push_back(vCount);           // Salva ID
                    	row.push_back(vCount); // Aggiungi indice del vertice alla riga corrente
                    	vCount++; // Avanza contatore vertice
					}
                }
                grid.push_back(row); // Aggiungi riga alla griglia
            }
            // variabile temporanea che memorizza gli estremi dei lati
            vector<Vector2i> eList;

            //Creiamo i nuovi lati dati dalla triangolazione e aggiorniamo la lista dei lati
            for(size_t i=0; i<grid.size(); i++){
                Vector2i extrema;
                for(size_t j=0; j<grid[i].size(); j++){
                    if(i<grid.size()-1){
                    	extrema << grid[i][j], grid[i + 1][j];
                    	
                    	if(!EdgeIsDupe(triMesh, extrema)){
                    		triMesh.Cell1DsId[eCount] = eCount;;
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  //lato sotto a sinistra
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i + 1][j];
	                        eList.push_back(extrema);
	                        eCount++;	
						}
						
						extrema << grid[i][j], grid[i + 1][j+1];
						if(!EdgeIsDupe(triMesh, extrema)){
	                        triMesh.Cell1DsId[eCount]=eCount;
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  //lato sotto a destra
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i + 1][j + 1];
	                        eList.push_back(extrema);
	                        eCount++;
	                    }
                    }
                    
                    if(i>0 && j<grid[i].size()-1){
                    	extrema << grid[i][j], grid[i][j+1];
                    	if(!EdgeIsDupe(triMesh, extrema)){
	                        triMesh.Cell1DsId[eCount] = eCount;;
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  //lato accanto a destra
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i][j + 1];
	                        eList.push_back(extrema);
	                        eCount++;
                    	}
                    }
                }
            }

			//Creiamo le nuove facce dopo la triangolazione
			for (size_t i=0; i<grid.size()-1; i++){
				for (size_t j=0; j<grid[i].size(); j++){
					triMesh.Cell2DsId[fCount]=(fCount);
					triMesh.Cell2DsVertices[fCount] = {grid[i][j],grid[i+1][j],grid[i+1][j+1]};
					for(unsigned int v=0; v<3; v++){
						unsigned int from=triMesh.Cell2DsVertices[fCount][v];
						unsigned int to=triMesh.Cell2DsVertices[fCount][(v+1)%3];
						vector<unsigned int> edges;
						for(unsigned int k=0; k<triDimensions[1]; k++){
							if((from==triMesh.Cell1DsExtrema(0, k) && to==triMesh.Cell1DsExtrema(1, k)) || (from==triMesh.Cell1DsExtrema(1, k) && to==triMesh.Cell1DsExtrema(0, k))){
								triMesh.Cell2DsEdges[fCount][v]=k;
								break;
							}
						}
					}
					fCount++;
					
					if(i>0 && j<grid[i].size()-1){
						triMesh.Cell2DsId[fCount]=fCount;
						triMesh.Cell2DsVertices[fCount] ={grid[i][j],grid[i][j+1],grid[i+1][j+1]};
						for(unsigned int v=0; v<3; v++){
							unsigned int from=triMesh.Cell2DsVertices[fCount][v];
							unsigned int to=triMesh.Cell2DsVertices[fCount][(v+1)%3];
							vector<unsigned int> edges;
							for(unsigned int k=0; k<triDimensions[1]; k++){
								if((from==triMesh.Cell1DsExtrema(0, k) && to==triMesh.Cell1DsExtrema(1, k)) || (from==triMesh.Cell1DsExtrema(1, k) && to==triMesh.Cell1DsExtrema(0, k))){
									triMesh.Cell2DsEdges[fCount][v]=k;
									break;
								}
							}
							
						}
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

		triMesh.Cell0DsCoordinates.conservativeResize(vCount, 3);
		return true;  
	}
	
	// Funzione per triangolare le facce e popolare la mesh triangolata di classe II
    bool GenerateTriangulatedMesh2(
        PolyhedralMesh& baseMesh,     // Mesh di partenza (con facce non triangolate)
        PolyhedralMesh& triMesh,      // Mesh di output triangolata
        const unsigned int& b, const unsigned int& c, // Parametri della suddivisione
        const Vector3i& triDimensions) // Dimensione di (V,E,F) della mesh triangolata (con duplicati)
    {
        unsigned int level = b + c; // Numero di suddivisioni laterali per triangolo

        // Inizializzazione della struttura dati della mesh triangolata

        // Allocazione memoria per vertici (0D)
        unsigned int guessVertices = 2 * triDimensions[0]; // Abbondante per duplicati
		triMesh.Cell0DsCoordinates = MatrixXd::Zero(guessVertices, 3); 
		triMesh.Cell0DsId.reserve(guessVertices);

        // Allocazione memoria per lati (1D)
        triMesh.Cell1DsId.resize(triDimensions[1]);
        triMesh.Cell1DsExtrema = MatrixXi::Zero(2, triDimensions[1]);

        // Allocazione memoria per facce (2D)
        triMesh.Cell2DsId.resize(triDimensions[2]);
        triMesh.Cell2DsEdges.resize(triDimensions[2]);
        for (auto& edgeList : triMesh.Cell2DsEdges) {
            edgeList.resize(3); // Ogni faccia ha 3 spigoli
        }
		for (auto& vertList : triMesh.Cell2DsVertices) {
        vertList.resize(3); // Ogni faccia ha 3 vertici
		}	
		
        // Id per vertici
        unsigned int vCount = 0;
		unsigned int eCount = 0;
		unsigned int fCount = 0;

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
					if(!VertexIsDupe(triMesh, pos)){
						triMesh.Cell0DsCoordinates.col(vCount) = pos; // Salva posizione
                    	triMesh.Cell0DsId.push_back(vCount);           // Salva ID
                    	row.push_back(vCount); // Aggiungi indice del vertice alla riga corrente
                    	vCount++; // Avanza contatore vertice
					}
                }
                grid.push_back(row); // Aggiungi riga alla griglia
            }
            // variabile temporanea che memorizza gli estremi dei lati
            vector<Vector2i> eList;

            //Creiamo i nuovi lati dati dalla triangolazione e aggiorniamo la lista dei lati
            for(size_t i=0; i<grid.size(); i++){
                Vector2i extrema;
                for(size_t j=0; j<grid[i].size(); j++){
                    if(i<grid.size()-1){
                    	extrema << grid[i][j], grid[i + 1][j+1];
                    	if(!EdgeIsDupe(triMesh, extrema)){
                    		triMesh.Cell1DsId[eCount] = eCount;;
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  //lato sotto a sinistra
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i + 1][j];
	                        eList.push_back(extrema);
	                        eCount++;	
						}
						
						extrema << grid[i][j], grid[i + 1][j+1];
						if(!EdgeIsDupe(triMesh, extrema)){
	                        triMesh.Cell1DsId[eCount]=eCount;
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  //lato sotto a destra
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i + 1][j + 1];
	                        eList.push_back(extrema);
	                        eCount++;
	                    }
                    }
                    
                    if(j<grid[i].size()-1){
                    	extrema << grid[i][j], grid[i + 1][j+1];
                    	if(!EdgeIsDupe(triMesh, extrema)){
	                        triMesh.Cell1DsId[eCount] = eCount;;
	                        triMesh.Cell1DsExtrema(0, eCount) = grid[i][j];  //lato accanto a destra
	                        triMesh.Cell1DsExtrema(1, eCount) = grid[i][j + 1];
	                        eList.push_back(extrema);
	                        eCount++;
                    	}
                    }
                }
            }
			
			//Creiamo le nuove facce dopo la triangolazione
			for (size_t i=0; i<grid.size(); i++){
				for (size_t j=0; j<grid[i].size(); j++){
					if (i<grid.size() -1) {
						triMesh.Cell2DsId[fCount]=(fCount);
						triMesh.Cell2DsVertices[fCount] = {grid[i][j],grid[i+1][j],grid[i+1][j+1]};
						for(unsigned int v=0; v<3; v++){
							unsigned int from=triMesh.Cell2DsVertices[fCount][v];
							unsigned int to=triMesh.Cell2DsVertices[fCount][(v+1)%3];
							vector<unsigned int> edges;
							for(Eigen::Index k=0; k<triMesh.Cell1DsExtrema.cols(); k++){
								if((from==triMesh.Cell1DsExtrema(0, k) && to==triMesh.Cell1DsExtrema(1, k)) || (from==triMesh.Cell1DsExtrema(1, k) && to==triMesh.Cell1DsExtrema(0, k))){
									triMesh.Cell2DsEdges[fCount][v]=k;
									break;
								}
							}
						}
						fCount ++;
					}
					
					if(j<grid[i].size()-1){
						triMesh.Cell2DsId[fCount]=fCount;
						triMesh.Cell2DsVertices[fCount] ={grid[i][j],grid[i][j+1],grid[i+1][j+1]};
						for(unsigned int v=0; v<3; v++){
							unsigned int from=triMesh.Cell2DsVertices[fCount][v];
							unsigned int to=triMesh.Cell2DsVertices[fCount][(v+1)%3];
							vector<unsigned int> edges;
							for(Eigen::Index k=0; k<triMesh.Cell1DsExtrema.cols(); k++){
								if((from==triMesh.Cell1DsExtrema(0, k) && to==triMesh.Cell1DsExtrema(1, k)) || (from==triMesh.Cell1DsExtrema(1, k) && to==triMesh.Cell1DsExtrema(0, k))){
									triMesh.Cell2DsEdges[fCount][v]=k;
									break;
								}
							}
							
						}
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

		triMesh.Cell0DsCoordinates.conservativeResize(vCount, 3);
		return true;  
	}	
}
