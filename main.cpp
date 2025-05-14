#include <iostream>
#include <Eigen/Eigen>
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

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
    vector<string> vector_input;
    string word;

    while (converter >> word) //finchè ci sono parole nello stream
    { 
        vector_input.push_back(word); //aggiungi la parola al vettore
    }
    
    if (vector_input.size() != 4 && vector_input.size() != 6) // controlla se il numero di elementi è 4 o 6
    {
        cerr << "Numero di elementi non valido!" << endl;
        return 1;
    }

    unsigned int p = vector_input[0];
    unsigned int q = vector_input[1];
    unsigned int b = vector_input[2];
    unsigned int c = vector_input[3];

    if (vector_input.size() == 6) // controlla se il numero di elementi è 6
    {
        unsigned int id_vertex_1 = vector_input[4];
        unsigned int id_vertex_2 = vector_input[5];
    }

    if(q!=3 && q!=4 && q!=5) // controlla se q è 3, 4 o 5
    {
        cerr << "q non valido!" << endl;
        return 1;
    }

    if(p!=3 && p!=4 && p!=5) // controlla se p è 3, 4 o 5
    {
        cerr << "p non valido!" << endl;
        return 1;
    }
    length_input = vector_input.size(); 


    // caso 1: q=3 p=3 export tetraedo
    // caso 2: q=3 p=4 export ottaedro
    // caso 3: q=3 p=5 export icosaedro
    // caso 4: q=4 p=3 export cubo (duale ottaedro)
    // caso 5: q=5 p=5 export dodecaedro (duale icosaedro)


}    