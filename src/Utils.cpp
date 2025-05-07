#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace PolygonalLibrary
{
bool ImportMesh(PolygonalMesh& mesh)
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
bool ImportCell0Ds(PolygonalMesh& mesh)
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
        replace(line.begin(), line.end(), ';', ' ');

        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2d coord;

        converter >>  id >> marker >> mesh.Cell0DsCoordinates(0, id) >> mesh.Cell0DsCoordinates(1, id);

        mesh.Cell0DsId.push_back(id);

        /// Memorizza i marker
        if(marker != 0)
        {
            const auto it = mesh.MarkerCell0Ds.find(marker);
            if(it == mesh.MarkerCell0Ds.end())
            {
                mesh.MarkerCell0Ds.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell0Ds[marker].push_back(id);
                it->second.push_back(id);
            }
        }

    }

    return true;
}
// ***************************************************************************
bool ImportCell1Ds(PolygonalMesh& mesh)
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
        replace(line.begin(), line.end(), ';', ' ');

        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2i extrema;

        converter >>  id >> marker >>  mesh.Cell1DsExtrema(0, id) >>  mesh.Cell1DsExtrema(1, id);
        mesh.Cell1DsId.push_back(id);

        /// Memorizza i marker
        if(marker != 0)
        {
            const auto it = mesh.MarkerCell1Ds.find(marker);
            if(it == mesh.MarkerCell1Ds.end())
            {
                mesh.MarkerCell1Ds.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell1Ds[marker].push_back(id);
                it->second.push_back(id);
            }
        }
    }

    return true;
}
// ***************************************************************************
bool ImportCell2Ds(PolygonalMesh& mesh)
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
        replace(line.begin(), line.end(), ';', ' ');

        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        unsigned int numVertices;
        unsigned int numEdges;

        converter >>  id >> marker >> numVertices;
        vector<unsigned int> vertices(numVertices);
        for(unsigned int i=0; i<numVertices; i++)
            converter >> vertices[i];
        converter >> numEdges;
        vector<unsigned int> edges(numEdges);
        for(unsigned int i = 0; i<numEdges; i++)
            converter >> edges[i];

        mesh.Cell2DsId.push_back(id);
        mesh.Cell2DsNumVertices.push_back(numVertices);
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsNumEdges.push_back(numEdges);
        mesh.Cell2DsEdges.push_back(edges);

        /// Memorizza i marker
        if(marker != 0)
        {
            const auto it = mesh.MarkerCell2Ds.find(marker);
            if(it == mesh.MarkerCell2Ds.end())
            {
                mesh.MarkerCell2Ds.insert({marker, {id}});
            }
            else
            {
                // mesh.MarkerCell2Ds[marker].push_back(id);
                it->second.push_back(id);
            }
        }
    }
    return true;
}

bool CheckMarkers(const PolygonalMesh& mesh)
{
    for (const auto& [marker, ids] : mesh.MarkerCell0Ds)
        for (const auto id : ids)
            if (std::find(mesh.Cell0DsId.begin(), mesh.Cell0DsId.end(), id) == mesh.Cell0DsId.end())
                return false;

    for (const auto& [marker, ids] : mesh.MarkerCell1Ds)
        for (const auto id : ids)
            if (std::find(mesh.Cell1DsId.begin(), mesh.Cell1DsId.end(), id) == mesh.Cell1DsId.end())
                return false;

    for (const auto& [marker, ids] : mesh.MarkerCell2Ds)
        for (const auto id : ids)
            if (std::find(mesh.Cell2DsId.begin(), mesh.Cell2DsId.end(), id) == mesh.Cell2DsId.end())
                return false;

    return true;
}

bool CheckEdges(const PolygonalMesh& mesh)
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

bool CheckAreas(const PolygonalMesh& mesh)
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
