#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
bool 
}




namespace PolyhedralLibrary
{
bool ImportMesh(PolyhedralMesh& mesh)
{

    if(!ImportCell0Ds(mesh))
        return false;

    if(!ImportCell1Ds(mesh))
        return false;

    if(!ImportCell2Ds(mesh))
        return false;

    return true;

}
// ***************************************************************************
bool ImportCell0Ds(PolyhedralMesh& mesh)
{
    ifstream file("./Cell0Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell0Ds = listLines.size();

    if (mesh.NumCell0Ds == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }

    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (string& line : listLines)
    {
        //replace(line.begin(), line.end(), ';', ' ');

        istringstream converter(line);

        unsigned int id;
        Vector3d coord;

        converter >>  id >> mesh.Cell0DsCoordinates(0, id) >> mesh.Cell0DsCoordinates(1, id) >> mesh.Cell0DsCoordinates(2, id);

        mesh.Cell0DsId.push_back(id);

    }

    return true;
}
// ***************************************************************************
bool ImportCell1Ds(PolyghedralMesh& mesh)
{
    ifstream file("./Cell1Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();

    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds);

    for (string& line : listLines)
    {
        //replace(line.begin(), line.end(), ';', ' ');

        istringstream converter(line);

        unsigned int id;
        Vector2i extrema;

        converter >>  id >>  mesh.Cell1DsExtrema(0, id) >>  mesh.Cell1DsExtrema(1, id);
        mesh.Cell1DsId.push_back(id);
    }

    return true;
}
// ***************************************************************************
bool ImportCell2Ds(PolyghedralMesh& mesh)
{
    ifstream file;
    file.open("./Cell2Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell2Ds = listLines.size();

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsNumVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsNumEdges.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    for (string& line : listLines)
    {
        //replace(line.begin(), line.end(), ';', ' ');

        istringstream converter(line);

        unsigned int id;
        unsigned int numVertices;
        unsigned int numEdges;

        converter >>  id >> numVertices >>  numEdges;
        vector<unsigned int> vertices(numVertices);
        for(unsigned int i=0; i<numVertices; i++)
            converter >> vertices[i];
        
        vector<unsigned int> edges(numEdges);
        for(unsigned int i = 0; i<numEdges; i++)
            converter >> edges[i];

        bool oriented = true;
        for (unsigned int e = 0; e < numEdges; e++) 
        {
            unsigned int currentEdge = edges[e];
            unsigned int nextEdge = edges[(e + 1) % numEdges];
        
            unsigned int nextEdgeOrigin = mesh.Cell1DsExtrema(0, nextEdge);
            unsigned int currentEdgeEnd = mesh.Cell1DsExtrema(1, currentEdge);
            unsigned int currentVertex = vertices[e];
            unsigned int currentEdgeOrigin = mesh.Cell1DsExtrema(0, currentEdge);

            // Controlla se il lato corrente finisce dove inizia il lato successivo
            if (currentEdgeEnd != nextEdgeOrigin || currentVertex != currentEdgeOrigin)
            {
                oriented = false;
                break;
            }
        }
        
        // Se l'orientamento non è corretto, prova a invertirlo
        if (!oriented) {
            reverse(vertices.begin(), vertices.end());
            reverse(edges.begin(), edges.end());
        }
        
        mesh.Cell2DsId.push_back(id);
        mesh.Cell2DsNumVertices.push_back(numVertices);
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsNumEdges.push_back(numEdges);
        mesh.Cell2DsEdges.push_back(edges);
    }

    return true;
}

bool ImportCell3Ds(PolyhedralMesh& mesh)
{
    ifstream file;
    file.open("./Cell3Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    // remove header
    listLines.pop_front();

    mesh.NumCell3Ds = listLines.size();

    if (mesh.NumCell3Ds == 0)
    {
        cerr << "There is no cell 3D" << endl;
        return false;
    }

    mesh.Cell3DsId.reserve(mesh.NumCell3Ds);
    mesh.Cell3DsNumVertices.reserve(mesh.NumCell3Ds);
    mesh.Cell3DsVertices.reserve(mesh.NumCell3Ds);
    mesh.Cell3DsNumEdges.reserve(mesh.NumCell3Ds);
    mesh.Cell3DsEdges.reserve(mesh.NumCell3Ds);
    mesh.Cell3DsNumFaces.reserve(mesh.NumCell3Ds);
    mesh.Cell3DsFaces.reserve(mesh.NumCell3Ds);

    for (string& line : listLines)
    {
        //replace(line.begin(), line.end(), ';', ' ');

        istringstream converter(line);

        unsigned int id;
        unsigned int numVertices;
        unsigned int numEdges;
        unsigned int numFaces;

        converter >>  id >> numVertices >>  numEdges >> numFaces;
        vector<unsigned int> vertices(numVertices);
        for(unsigned int i=0; i<numVertices; i++)
            converter >> vertices[i];
        vector<unsigned int> edges(numEdges);
        for(unsigned int i = 0; i<numEdges; i++)
            converter >> edges[i];
        vector<unsigned int> faces(numFaces);
        for(unsigned int i = 0; i<numFaces; i++)
            converter >> faces[i];

        mesh.Cell3DsId.push_back(id);
        mesh.Cell3DsNumVertices.push_back(numVertices);
        mesh.Cell3DsVertices.push_back(vertices);
        mesh.Cell3DsNumEdges.push_back(numEdges);
        mesh.Cell3DsEdges.push_back(edges);
        mesh.Cell3DsNumFaces.push_back(numFaces);
        mesh.Cell3DsFaces.push_back(faces);
    }
    return true;
}

bool CheckEdges(const PolyghedralMesh& mesh)
{
    for (unsigned int i = 0; i < mesh.NumCell1Ds; i++)
    {
        unsigned int id = mesh.Cell1DsId[i];
        unsigned int origin = mesh.Cell1DsExtrema(0, id);
        unsigned int end = mesh.Cell1DsExtrema(1, id);

        Vector2d A = mesh.Cell0DsCoordinates.col(origin);
        Vector2d B = mesh.Cell0DsCoordinates.col(end);
        double length = (A - B).norm();

        if (length == 0)
        {
            cerr << "Null edge" << i << endl;
            return false;
        }
    }
    return true;
}

bool CheckAreas(const PolyghedralMesh& mesh)
{
    for (unsigned int i = 0; i < mesh.NumCell2Ds; i++)
    {
        const auto& vertices = mesh.Cell2DsVertices[i];
        double area = 0.0;

        for (unsigned int j = 0; j < vertices.size(); j++)
        {
            Vector2d p1 = mesh.Cell0DsCoordinates.col(vertices[j]);
            Vector2d p2 = mesh.Cell0DsCoordinates.col(vertices[(j + 1) % vertices.size()]);
            area += p1(0) * p2(1) - p2(0) * p1(1);
        }

        area = abs(area) / 2.0;

        if (area == 0)
        {
            cerr << "Null area " << i << endl;
            return false;
        }
    }

    return true;
}

}
