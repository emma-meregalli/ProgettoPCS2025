 #pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolyhedralLibrary {

struct PolyhedralMesh
{
    // Vertici
    unsigned int NumCell0Ds = 0;    // number of Cell0Ds
    std::vector<unsigned int> Cell0DsId = {};   // Cell0Ds id
    Eigen::MatrixXd Cell0DsCoordinates = {};    // Cell0Ds coordinates
    vector<bool> Cell0DsDupes = {}; // Cell0Ds duplicati

    // Spigoli
    unsigned int NumCell1Ds = 0;    // number of Cell1Ds
    std::vector<unsigned int> Cell1DsId = {};   // Cell1Ds id
    Eigen::MatrixXi Cell1DsExtrema = {};    // Cell1Ds vertices id (start, end)
    vector<bool> Cell1DsDupes = {}; // Cell1Ds duplicati

    //Facce
    unsigned int NumCell2Ds = 0;    // number of Cell2D
    std::vector<unsigned int> Cell2DsId = {};   // Cell2D id
	std::vector<unsigned int> Cell2DsNumVertices = {};
	std::vector<unsigned int> Cell2DsNumEdges = {};
    std::vector<std::vector<unsigned int>> Cell2DsVertices; // Cell2D vertices id
    std::vector<std::vector<unsigned int>> Cell2DsEdges = {};   // Cell2D edges id

    //Poliedri
    unsigned int NumCell3Ds = 0;    // number of Cell3D
    std::vector<unsigned int> Cell3DsId = {};   // Cell3D id
    std::vector<unsigned int> Cell3DsVertices = {};     // Cell3D vertices id
    std::vector<unsigned int> Cell3DsEdges = {};   // Cell3D edges id
    std::vector<unsigned int> Cell3DsFaces = {};    // Cell3D faces id
	
	//Cammini minimi
	std::vector<unsigned int> Cell0DsShortPath;
	std::vector<unsigned int> Cell1DsShortPath;
	std::vector<unsigned int> Cell1DsEsistente;
};

}
