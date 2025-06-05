#include <iostream>
#include <Eigen/Eigen>
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int checkInput(const vector<int> vector_input){
    unsigned int p = vector_input[0];
    unsigned int q = vector_input[1];
    unsigned int b = vector_input[2];
    unsigned int c = vector_input[3];

    PolyhedralMesh mesh;
    PolyhedralMesh triMesh;
    PolyhedralMesh dualMesh;
    

	if(p!=3 && p!=4 && p!=5) // controlla se p è 3, 4 o 5
    {
        cerr << "p non valido!" << endl;
        return 1;
    }
    if(q!=3 && q!=4 && q!=5) // controlla se q è 3, 4 o 5
    {
        cerr << "q non valido!" << endl;
        return 1;
    }
    
    if(b!=c && b!=0 && c!=0){
    	cerr << "b e c non validi!" <<endl;
	}

    if (b==0 || c==0){
        if (p==3){
            switch (q){
                case 3 : //genera il tetraedro se p=3,q=3
                    ExportTetrahedron(mesh,triMesh,b,c);
                    CreateTxtFiles(triMesh);
                case 4 : //genera l'ottaedro se p=3,q=4
                    ExportOctahedron(mesh,triMesh,b,c);
                    CreateTxtFiles(triMesh);
                case 5 : //genera l'icosaedro se p=3,q=5
                    ExportIcosahedron(mesh,triMesh,b,c);
                    CreateTxtFiles(triMesh);
            }
        if (q==3){
            switch (p){
                case 3 : //genera il tetraedro se p=3,q=3
                    ExportTetrahedron(mesh,triMesh,b,c);
                    CreateTxtFiles(triMesh);
                case 4 : //genera il tetraedro se p=4,q=3
                    ExportOctahedron(mesh,triMesh,b,c);
                    GenerateDual(triMesh, dualMesh);
                    CreateTxtFiles(dualMesh);
                case 5 : //genera il tetraedro se p=5,q=3
                    ExportIcosahedron(mesh,triMesh,b,c);
                    GenerateDual(triMesh, dualMesh);
                    CreateTxtFiles(dualMesh);
            }
        }
        }
    }
    
    if(b==c){
    	//Poliedri geodetici di classe II
	}
	return 1
}

int main()
{
    cout << "Inserire una quadrupla di numeri interi del tipo (p,q,b,0) oppure (p,q,0,c), altrimenti inserire una sestupla del tipo (p,q,b,c, id_vertice_1, id_vertice_2)"<< endl;
    string input = "";
    cin >> input;

    if (input.front() == '(' && input.back() == ')') // controlla se il formato dell'input è corretto
    {
        input = input.substr(1, input.size() - 2); //crea una sottostringa a partire dal carattere 1(che sarebbe quello dopo "(" di lunghezza la lunghezza della stringa -2
    } 
    else 
    {
        cerr << "Formato non valido!" << endl;
        return 1;
    }

    replace(input.begin(), input.end(), ',', ' ');
    istringstream converter(input);  // crea uno stream di input dalla stringa
    vector<int> vector_input;
    string word;

    while (converter >> word) //finchè ci sono parole nello stream
    { 
        vector_input.push_back(stoi(word)); //aggiungi la parola al vettore
    }
    
    if (vector_input.size() != 4 && vector_input.size() != 6) // controlla se il numero di elementi è 4 o 6
    {
        cerr << "Numero di elementi non valido!" << endl;
        return 1;
    }

    if (vector_input.size() == 4)
        checkInput(vector_input);

    if (vector_input.size() == 6) // controlla se il numero di elementi è 6
    {
        unsigned int id_vertex_2 = vector_input[5];
        vector_input.pop_back();
        unsigned int id_vertex_1 = vector_input[4];
        vector_input.pop_back();

        checkInput(vector_input);
    }
}    
