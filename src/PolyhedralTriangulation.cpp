#include <vector>  
#include <array>  
#include <Eigen/Dense>  
#include <limits>       
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
    
    //AddOrFindEdge ha lo scopo di evitare duplicati durante la costruzione della mesh triangolata, aggiungendo un lato (edge) solo se non esiste già 

    // Funzione ausiliaria per aggiungere un lato alla mesh triangolata, evitando duplicati
    void AddOrFindEdge(
        int vStart, int vEnd, // estremi del lato
        PolyhedralMesh& triMesh, 
        unsigned int& edgeCounter, // contatore globale dei lati
        unsigned int faceIndex) // indice della faccia corrente nella mesh
    {
        bool alreadyExists = false; // variabile booleana per indicare se il lato esiste già

        // Verifico se il lato da vStart a vEnd esiste già nella mesh
        for (unsigned int i = 0; i < edgeCounter; i++)  //Scorro tutti i lati già inseriti nella mesh triangolata
        {
            if ((triMesh.Cell1DsExtrema(i, 0) == vStart && triMesh.Cell1DsExtrema(i, 1) == vEnd) ||  
                (triMesh.Cell1DsExtrema(i, 0) == vEnd && triMesh.Cell1DsExtrema(i, 1) == vStart)) //Controllo se l’i-esimo lato ha come estremi vStart e vEnd, nell'ordine dato e nell'ordine inverso  
            {
                triMesh.Cell2DsEdges[faceIndex].push_back(i); //Se il lato esiste già, viene aggiunto alla lista dei lati della faccia corrente (faceIndex).
                alreadyExists = true; //setto alreadyExist = true per evitare di aggiungere nuovamente il lato.
                break;
            }
        }

        if (!alreadyExists) // se il lato non esiste già lo aggiungo alla mesh triangolata
        {
            triMesh.Cell1DsExtrema.row(edgeCounter) << vStart, vEnd; // Inserisco gli estremi del lato nella tabella dei lati 
            triMesh.Cell1DsId[edgeCounter] = edgeCounter; // Assegno l’ID del lato (uguale all’indice in cui viene inserito)

            const auto& flagsA = triMesh.Cell0DsFlag[vStart]; // Flag associato al primo vertice (vStart)
            const auto& flagsB = triMesh.Cell0DsFlag[vEnd];   // Flag associato al secondo vertice (vEnd)
            //I flag dicono a quali facce appartiene ciascun vertice

            bool hasCommonFlag = false; // variabile booleana per verificare se i due vertici hanno flag comuni
            // Cerco i flag comuni tra i due vertici
            for (unsigned int fA : flagsA) // Scorro i flag del primo vertice
            {
                for (unsigned int fB : flagsB)  // Scorro i flag del secondo vertice
                {
                    if (fA == fB) // Se trovo un flag comune tra i due vertici
                    {
                        triMesh.Cell1DsFlag[edgeCounter] = fA; // Assegno il flag comune al lato
                        // Se i flag sono uguali, significa che il lato appartiene a quella faccia
                        hasCommonFlag = true;
                        break;
                    }
                }
                if (hasCommonFlag) break; // Se ho trovato un flag comune, esco dal ciclo
            }

            if (!hasCommonFlag) // Se non ci sono flag comuni, assegno un valore di default invalido
            {
                triMesh.Cell1DsFlag[edgeCounter] = numeric_limits<unsigned int>::max(); // Se non hanno nulla in comune, assegno flag "invalido"
            }
            triMesh.Cell2DsEdges[faceIndex].push_back(edgeCounter); // Aggiungo l'indice del lato alla lista dei lati della faccia corrente
            edgeCounter++; // Incremento contatore globale dei lati
        }
    }

    // Funzione principale per triangolare le facce e popolare la mesh triangolata (SONO QUIIIIIIIIIIII)
    void GenerateTriangulatedMesh(
        PolyhedralMesh& baseMesh,     // Mesh di partenza (con facce non triangolate)
        PolyhedralMesh& triMesh,      // Mesh di output triangolata
        unsigned int b, unsigned int c, // Parametri della suddivisione
        const vector<int>& duplicatedDimensions) // Dimensione di (V,E,F) della mesh triangolata (con duplicati)
        {
        unsigned int refLevel = b + c; // Numero di suddivisioni laterali per triangolo

        // Inizializzazione della struttura dati della mesh triangolata

        // Allocazione memoria per vertici (0D)
        triMesh.Cell0DsId.resize(duplicatedDimensions[0]);
        triMesh.Cell0DsCoordinates = MatrixXd::Zero(duplicatedDimensions[0],3);
        triMesh.Cell0DsFlag.resize(duplicatedDimensions[0]);

        // Allocazione memoria per lati (1D)
        triMesh.Cell1DsId.resize(duplicatedDimensions[1]);
        triMesh.Cell1DsExtrema = MatrixXi::Zero(duplicatedDimensions[1], 2);
        triMesh.Cell1DsFlag.resize(duplicatedDimensions[1]);

        // Allocazione memoria per facce (2D)
        triMesh.Cell2DsId.resize(duplicatedDimensions[2]);
        triMesh.Cell2DsVertices.resize(duplicatedDimensions[2],3);
        triMesh.Cell2DsEdges.resize(duplicatedDimensions[2],3);

        // Contatori per vertici, lati e facce
        unsigned int vCount = 0, eCount = 0, fCount = 0;

        // Ciclo su tutte le facce della mesh di base
        for (unsigned int faceIdx = 0; faceIdx < baseMesh.Cell2DsId.size(); faceIdx++) 
        {
            const auto& faceVerts = baseMesh.Cell2DsVertices[faceIdx]; //Prendo i tre vertici della faccia corrente
        
            // Coordinate dei 3 vertici del triangolo originale
            Vector3d A = baseMesh.Cell0DsCoordinates.col(faceVerts[0]); // Vertice A
            Vector3d B = baseMesh.Cell0DsCoordinates.col(faceVerts[1]); // Vertice B
            Vector3d C = baseMesh.Cell0DsCoordinates.col(faceVerts[2]); // Vertice C

            vector<vector<int>> grid; // Griglia di vertici interni alla faccia
            // la griglia ha refLevel = b + c righe e ogni riga i ha i + 1 elementi (forma triangolare)
            // Costruzione della griglia interplata sulla faccia
            for (unsigned int i = 0; i <= refLevel; i++) {
                vector<int> row; //riga corrente dei vertici 

                // Calcolo il punto iniziale e finale della riga i-esima
                Vector3d from = ((double)i / refLevel) * B + ((double)(refLevel - i) / refLevel) * A;
                Vector3d to = ((double)i / refLevel) * C + ((double)(refLevel - i) / refLevel) * A;

                for (unsigned int j = 0; j <= i; j++) {
                    // Interpola tra from e to per ottenere un punto interno
                    Vector3d pos;
                    if (i == 0) {
                        pos = A;
                    } else {
                        pos = ((double)j / i) * to + ((double)(i - j) / i) * from;
                        pos.normalize(); // Proietta sulla sfera unitaria
                    }
                    triMesh.Cell0DsCoordinates.col(vCount) = pos; // Salva posizione
                    triMesh.Cell0DsId[vCount] = vCount;           // Salva ID

                    // Assegna flag in base alla posizione nella griglia
                    if (i == 0) 
                    {
                        triMesh.Cell0DsFlag[vCount] = {baseMesh.Cell2DsEdges[faceIdx][0], baseMesh.Cell2DsEdges[faceIdx][2]};
                    } else if (i == refLevel) {
                        if (j == 0)
                            triMesh.Cell0DsFlag[vCount] = {baseMesh.Cell2DsEdges[faceIdx][0], baseMesh.Cell2DsEdges[faceIdx][1]};
                        else if (j == refLevel)
                            triMesh.Cell0DsFlag[vCount] = {baseMesh.Cell2DsEdges[faceIdx][1], baseMesh.Cell2DsEdges[faceIdx][2]};
                        else
                            triMesh.Cell0DsFlag[vCount] = {baseMesh.Cell2DsEdges[faceIdx][1]};
                    } else if (j == 0) {
                        triMesh.Cell0DsFlag[vCount] = {baseMesh.Cell2DsEdges[faceIdx][0]};
                    } else if (j == i) {
                        triMesh.Cell0DsFlag[vCount] = {baseMesh.Cell2DsEdges[faceIdx][2]};
                    } else {
                        triMesh.Cell0DsFlag[vCount] = {numeric_limits<unsigned int>::max()};
                    }

                    row.push_back(vCount); // Aggiungi indice del vertice alla riga corrente
                    vCount++; // Avanza contatore vertici
                }
                grid.push_back(row); // Aggiungi riga alla griglia
            }

            // Costruzione dei triangoli nella griglia (2 triangoli per quadrilatero)
            for (unsigned int i = 0; i < refLevel; i++) {
                for (unsigned int j = 0; j < i; j++) {
                    // Primo triangolo
                    vector<unsigned int> tri1 = {grid[i][j], grid[i + 1][j], grid[i + 1][j + 1]};
                    triMesh.Cell2DsVertices[fCount] = tri1;
                    triMesh.Cell2DsId[fCount] = fCount;
                    for (unsigned int e = 0; e < 3; e++)
                        AddOrFindEdge(tri1[e], tri1[(e + 1) % 3], triMesh, eCount, fCount);
                    fCount++;

                    // Secondo triangolo
                    vector<unsigned int> tri2 = {grid[i][j], grid[i + 1][j + 1], grid[i][j + 1]};
                    triMesh.Cell2DsVertices[fCount] = tri2;
                    triMesh.Cell2DsId[fCount] = fCount;
                    for (unsigned int e = 0; e < 3; e++)
                        AddOrFindEdge(tri2[e], tri2[(e + 1) % 3], triMesh, eCount, fCount);
                    fCount++;
                }

                // Triangolo finale sul bordo destro
                vector<unsigned int> lastTri = {grid[i][i], grid[i + 1][i], grid[i + 1][i + 1]};
                triMesh.Cell2DsVertices[fCount] = lastTri;
                triMesh.Cell2DsId[fCount] = fCount;
                for (unsigned int e = 0; e < 3; e++)
                    AddOrFindEdge(lastTri[e], lastTri[(e + 1) % 3], triMesh, eCount, fCount);
                fCount++;
            }
        }
    }

}
