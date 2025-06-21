 #pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolyhedralLibrary {

struct PolyhedralMesh
{
    // Vertici
    unsigned int NumCell0Ds = 0;    // Numero di elementi di Cell0Ds
    std::vector<unsigned int> Cell0DsId = {};   // Cell0Ds id
    Eigen::MatrixXd Cell0DsCoordinates = {};    // Cell0Ds coordinate

    // Spigoli
    unsigned int NumCell1Ds = 0;    // Numero di elementi di Cell1Ds
    std::vector<unsigned int> Cell1DsId = {};   // Cell1Ds id
    Eigen::MatrixXi Cell1DsExtrema = {};    // Cell1Ds id dei vertici (start, end)

    //Facce
    unsigned int NumCell2Ds = 0;    // Numero di elementi di Cell2D
    std::vector<unsigned int> Cell2DsId = {};   // Cell2D id
	std::vector<unsigned int> Cell2DsNumVertices = {};    // Numero di vertici della faccia
	std::vector<unsigned int> Cell2DsNumEdges = {};    // Numero di spigoli della faccia
    std::vector<std::vector<unsigned int>> Cell2DsVertices;    // Cell2D id dei vertici
    std::vector<std::vector<unsigned int>> Cell2DsEdges = {};    // Cell2D id degli spigoli

    //Poliedri
    unsigned int NumCell3Ds = 0;    // Numero di elementi di Cell3D
    std::vector<unsigned int> Cell3DsId = {};   // Cell3D id
    std::vector<unsigned int> Cell3DsVertices = {};     // Cell3D id dei vertici
    std::vector<unsigned int> Cell3DsEdges = {};   // Cell3D id degli spigoli
    std::vector<unsigned int> Cell3DsFaces = {};    // Cell3D id delle facce
	
	//Cammini minimi
	std::vector<unsigned int> Cell0DsShortPath;    // Lista dei vertici del cammino minimo
	std::vector<unsigned int> Cell1DsShortPath;    // Lista degli spigoli del cammino minimo
};

}
