#include <vector>  
#include <array>  
#include <Eigen/Dense>  
#include <limits> 
#include <utility>       
#include "Utils.hpp"   
#include "PolyhedralMesh.hpp" 

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {

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

    bool VertexIsDupe(const PolygonalMesh& mesh, const Vector3d& v){
        //Fisso una tolleranza per confrontare i vertici
	    double tol=1e-12;

        // Confronto con tutti i vertici già inseriti nella lista
        for (size_t i = 0; i < mesh.Cell0DsId.size(); i++) {
            if ((mesh.Cell0DsCoordinates.row(i) - v).norm() < tol) {    //Se il vertice esiste, allora restituisco il suo ID (bisogna fare un controllo con la tolleranza?)
                return true;  
            }
        }
        return false;
    }

    bool EdgeIsDupe(const PolygonalMesh& mesh, const Vector2i& e){
        for(size_t i=0; i<mesh.Cell1DsId.size(); i++){
            if(mesh.Cell1DsExtrema[i]==e)
                return true;
            swap(e[0],e[1]);
            if(mesh.Cell1DsExtrema[i]==e)
                return true;
        }
        return false;
    }
    
    // Funzione principale per triangolare le facce e popolare la mesh triangolata
    void GenerateTriangulatedMesh(
        PolyhedralMesh& baseMesh,     // Mesh di partenza (con facce non triangolate)
        PolyhedralMesh& triMesh,      // Mesh di output triangolata
        unsigned int b, unsigned int c, // Parametri della suddivisione
        const vector<int>& triDimensions) // Dimensione di (V,E,F) della mesh triangolata (con duplicati)
    {
        unsigned int level = b + c; // Numero di suddivisioni laterali per triangolo

        // Inizializzazione della struttura dati della mesh triangolata

        // Allocazione memoria per vertici (0D)
        triMesh.Cell0DsId.resize(triDimensions[0]);
        triMesh.Cell0DsCoordinates = MatrixXd::Zero(triDimensions[0],3);
        triMesh.Cell0DsFlag.resize(triDimensions[0]);

        // Allocazione memoria per lati (1D)
        triMesh.Cell1DsId.resize(triDimensions[1]);
        triMesh.Cell1DsExtrema = MatrixXi::Zero(triDimensions[1], 2);
        triMesh.Cell1DsFlag.resize(triDimensions[1]);

        // Allocazione memoria per facce (2D)
        triMesh.Cell2DsId.resize(triDimensions[2]);
        triMesh.Cell2DsVertices.resize(triDimensions[2],3);
        triMesh.Cell2DsEdges.resize(triDimensions[2],3);

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

            vector<vector<int>> grid; // Griglia di vertici interni alla faccia
            // la griglia ha level = b + c righe e ogni riga i ha i + 1 elementi (forma triangolare)
            // Costruzione della griglia interplata sulla faccia
            for (unsigned int i = 0; i <= level; i++) {
                vector<int> row; //riga corrente dei vertici 

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
                    triMesh.Cell0DsCoordinates.col(vCount) = pos; // Salva posizione
                    triMesh.Cell0DsId[vCount] = vCount;           // Salva ID

                    triMesh.Cell0DsDupes[vCount]=VertexIsDupe(triMesh, pos);  //restituisce True se il vertice esiste già nella lista

                    row.push_back(vCount); // Aggiungi indice del vertice alla riga corrente
                    vCount++; // Avanza contatore vertici
                }
                grid.push_back(row); // Aggiungi riga alla griglia
            }

            // variabile temporanea che memorizza gli estremi dei lati
            vector<Vector2i> eList;

            //Creiamo i nuovi lati dati dalla triangolazione e aggiorniamo la lista dei lati
            for(size_t i=0; i<grid.size(); i++){
                Vector2i estrema;
                for(size_t j=0; j<grid[i].size(); j++){
                    if(i<grid.size()-1){
                        triMesh.Cell1DsId[eCount]=eCount;
                        triMesh.Cell1DsExtrema[eCount]=[grid[i][j],grid[i+1][j]];  //lato sotto a sinistra
                        extrema=triMesh.Cell1DsExtrema[eCount];
                        triMesh.Cell1DsDupes[eCount]=EdgeIsDupe(triMesh, extrema);
                        eList.push_back(estrema);
                        eCount++;

                        triMesh.Cell1DsId[eCount]=eCount;
                        triMesh.Cell1DsExtrema[eCount]=[grid[i][j],grid[i+1][j+1]];  // lato sotto a destra
                        extrema=triMesh.Cell1DsExtrema[eCount];
                        triMesh.Cell1DsDupes[eCount]=EdgeIsDupe(triMesh, extrema);
                        eList.push_back(estrema);
                        eCount++;
                    }
                    if(j<grid[i].size()-1){
                        triMesh.Cell1DsId[eCount]=eCount;
                        triMesh.Cell1DsExtrema[eCount]=[grid[i][j],grid[i][j+1]];
                        extrema=triMesh.Cell1DsExtrema[eCount];
                        triMesh.Cell1DsDupes[eCount]=EdgeIsDupe(triMesh, extrema);
                        eList.push_back(estrema);
                        eCount++;
                    }
                }
            }
			
			unsigned int fCount = 0;
			for (size_t i=0; i<grid.size(); i++){
				for (size_t j=0; j<grid[i].size(); j++){
					if (i<grid.size() -1) {
						triMesh.Cell2DsId[fCount] = fCount;
						triMesh.Cell2DsVertices=[grid[i][j],grif[i+1][j],grid[i+1][j+1]];
						triMesh.Cell2DsEdges=[triMesh.Cell1DsId[  //dobbiamo aggiungere la memorizzazione dei lati della triangolazione che formano le facce
						
						
					}
					
						
			    }
		    }
		}

    //"nome funzione" ha lo scopo di popolare la cell3Ds
    
    }
}
