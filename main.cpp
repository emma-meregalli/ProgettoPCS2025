#include <iostream>
#include <Eigen/Eigen>
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int checkInput(const vector<int> vector_input){
	// Salva l'input dell'utente 
    unsigned int p = vector_input[0];
    unsigned int q = vector_input[1];
    unsigned int b = vector_input[2];
    unsigned int c = vector_input[3];
	
	// Per il cammino minimo, se il poliedro non viene triangolato (quindi se b+c==1) allora devo considerare i lati iniziali.
	// In caso contrario, devo considerare tutti i lati che si formano dopo la triangolazione
	
	//Controlla se serve fare la triangolazione
	bool all_edges = false;
	if (b + c == 1){
		all_edges = true;
	}	
	
	// Controlla se l'input prevede la ricerca del cammino minimo
	bool path = false;
	
	unsigned int id_v1 = 0;
	unsigned int id_v2 = 0;
	if (vector_input.size() == 6){
		id_v1 = vector_input[4];
		id_v2 = vector_input[5];
		path = true;
	}	

	// Mesh da usare nel programma
    PolyhedralMesh mesh;
    PolyhedralMesh triMesh;
    PolyhedralMesh dualMesh;
	bool shortestPath;
    
	// Controlla se p è uguale 3, 4 o 5
	if(p != 3 && p != 4 && p != 5){
        cerr << "p non valido!" << endl;
        return 1;
    }
    
    // Controlla se q è uguale 3, 4 o 5
    if(q != 3 && q != 4 && q != 5){
        cerr << "q non valido!" << endl;
        return 1;
    }
    
    // Controlla se b!=c e b,c!=0, che sono input non validi
    if(b != c && b != 0 && c != 0){
    	cerr << "b e c non validi!" <<endl;
    	return 1;
	}
	
	// Controlla se b=c=0, che sono input non validi
	if(b == 0 && c == 0){
    	cerr << "b e c non validi!" <<endl;
    	return 1;
	}

	// Decide quale poliedro va generato in base all'input
    if (p == 3){
        switch (q){
           case 3 : // Genera il tetraedro se p=3,q=3
                ExportTetrahedron(mesh, triMesh, b, c);
                cout << "Triangolato!" << endl;
                CreateTxtFiles(triMesh);
                cout << "Stampato!" << endl;
				if (path){
					shortestPath = ShortestPath(triMesh, id_v1, id_v2, triMesh.NumCell1Ds, all_edges);
					cout << "Cammino minimo trovato!" <<endl;
				}	
				ExportParaView(triMesh, path);
                break;
                
            case 4 : // Genera l'ottaedro se p=3,q=4
                ExportOctahedron(mesh, triMesh, b, c);
				cout << "Triangolato!" << endl;
                CreateTxtFiles(triMesh);
				cout << "Stampato!" << endl;
				if (path){
					shortestPath = ShortestPath(triMesh, id_v1, id_v2, triMesh.NumCell1Ds, all_edges);
					cout << "Cammino minimo trovato!" <<endl;
				}	
				ExportParaView(triMesh, path);
                break;
                
            case 5 : //genera l'icosaedro se p=3,q=5
                ExportIcosahedron(mesh, triMesh, b, c);
				cout << "Triangolato!" << endl;
                CreateTxtFiles(triMesh);
				cout << "Stampato!" << endl;
				if (path){
					shortestPath = ShortestPath(triMesh, id_v1, id_v2, triMesh.NumCell1Ds ,all_edges);
					cout << "Cammino minimo trovato!" << endl;
				}	
				ExportParaView(triMesh, path);
                break;
        }
    }
    if (q == 3 && p != 3){
        switch (p){
            case 4 : // Genera il tetraedro se p=4,q=3
                ExportOctahedron(mesh, triMesh, b, c);
				cout << "Triangolato!" << endl;
                GenerateDual(triMesh, dualMesh);
				cout << "Fatto il duale!" << endl;
                CreateTxtFiles(dualMesh);
				cout << "Stampato!" << endl;
				if (path){
					shortestPath = ShortestPath(dualMesh, id_v1, id_v2, dualMesh.NumCell1Ds ,all_edges);
					cout << "Cammino minimo trovato!" << endl;
				}	
				ExportParaView(dualMesh, path);
                break;
                
            case 5 : // Genera il tetraedro se p=5,q=3
                ExportIcosahedron(mesh, triMesh, b, c);
				cout << "Triangolato!" << endl;
                GenerateDual(triMesh, dualMesh);
				cout << "Fatto il duale!" << endl;
                CreateTxtFiles(dualMesh);
				cout << "Stampato!" << endl;
				if (path){
					shortestPath = ShortestPath(dualMesh, id_v1, id_v2, dualMesh.NumCell1Ds,all_edges);
					cout << "Cammino minimo trovato!" << endl;
				}	
				ExportParaView(dualMesh, path);
                break;
        }
    }
	return 0;	
}

int main(){
    cout << "Inserire una quadrupla di numeri interi del tipo (p,q,b,0), (p,q,0,c) oppure (p,q,b,b). Per i cammini minimi, inserire una sestupla del tipo (p,q,b,c,id_vertice_1,id_vertice_2) che segua le indicazioni precedenti."<< endl;
    string input = "";
    cin >> input;

	// Controlla se il formato dell'input è corretto
    if (input.front() == '(' && input.back() == ')'){
        input = input.substr(1, input.size() - 2); // Crea una sottostringa a partire dal carattere in posizione 1 (che sarebbe quello dopo "(" ) fino al penultimo carattere in posizione input.size()-2
    } 
    else{
        cerr << "Formato non valido!" << endl;
        return 1;
    }

    replace(input.begin(), input.end(), ',', ' '); //Toglie le parentesi, le virgole e gli spazi dalla stringa
    istringstream converter(input);  // Crea uno stream di input dalla stringa
    vector<int> vector_input;
    string word;

    while (converter >> word){ // Finchè ci sono parole nello stream 
        vector_input.push_back(stoi(word)); // Aggiunge la parola al vettore
    }
    
    // Controlla se il numero di elementi è 4 o 6
    if (vector_input.size() != 4 && vector_input.size() != 6){ 
        cerr << "Numero di elementi non valido!" << endl;
        return 1;
    }

	checkInput(vector_input);
    return 0;
}    
