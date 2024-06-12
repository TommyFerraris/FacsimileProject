#include "Triangolazione.hpp"

#include <fstream>

namespace GeometryLibrary{
//*********************************************************
vector<vector<vector<unsigned int>>> Polygons::TriangulatePolygons()
{
    const unsigned int numPolygons = listVertices.size();
    vector<vector<vector<unsigned int>>> triangleList(numPolygons);

    for(unsigned int p = 0; p < numPolygons; p++)
    {
        const unsigned int numPolygonVertices = listVertices[p].size();

        for (unsigned int v = 0; v < numPolygonVertices; v++)
        {
            const unsigned int nextVertex = listVertices[p][(v + 1) % numPolygonVertices];
            const unsigned int nextNextVertex = listVertices[p][(v + 2) % numPolygonVertices];

            if ((v + 2) % numPolygonVertices == 0)
                break;

            vector<unsigned int> triangle_vertices = {listVertices[p][0], nextVertex, nextNextVertex};

            triangleList[p].push_back(triangle_vertices);
        }
    }
    return triangleList;
}

void Polygons::GedimInterface(vector<vector<unsigned int>>& triangles,
                              VectorXi& materials)
{
    const unsigned int numPolygons = listVertices.size();
    vector<vector<vector<unsigned int>>> triangleList = TriangulatePolygons();

    unsigned int numTotalTriangles = 0;
    for(unsigned int p = 0; p < numPolygons; p++)
        numTotalTriangles += triangleList[p].size();

    triangles.reserve(numTotalTriangles);
    materials = VectorXi::Zero(numTotalTriangles);

    unsigned int count = 0;
    for(unsigned int p = 0; p < numPolygons; p++)
    {
        for(unsigned int t = 0; t < triangleList[p].size(); t++)
        {
            triangles.push_back(triangleList[p][t]);
            materials(count) = p;
            count++;
        }
    }
}
}
