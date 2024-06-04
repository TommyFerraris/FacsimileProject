#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace DFN_PolygonalLibrary {

struct PolygonalMesh
{
    unsigned int NumberCell0D = 0; ///< number of Cell0D
    vector<unsigned int> Cell0DId = {}; ///< Cell0D id, size 1 x NumberCell0D
    vector<Vector3d> Cell0DCoordinates = {}; ///< Cell0D coordinates, size 2 x NumberCell0D (x,y)

    unsigned int NumberCell1D = 0; ///< number of Cell1D
    vector<unsigned int> Cell1DId = {}; ///< Cell1D id, size 1 x NumberCell1D
    vector<Vector2i> Cell1DVertices = {}; ///< Cell1D vertices indices, size 2 x NumberCell1D (fromId,toId)

    unsigned int NumberCell2D = 0; ///< number of Cell2D
    vector<unsigned int> Cell2DId = {}; ///< Cell2D id, size 1 x NumberCell2D
    vector<unsigned int> Cell2DNumVertices = {}; ///< Cell2D number of vertices, size 1 x NumberCell2D
    vector<vector<unsigned int>> Cell2DVertices = {}; ///< Cell2D Vertices indices, size NumberCell2DVertices x NumberCell2D
    vector<unsigned int> Cell2DNumEdges = {}; ///< Cell2D number of vertices, size 1 x NumberCell2D
    vector<vector<unsigned int>> Cell2DEdges = {}; ///< Cell2D Cell1D indices, size NumberCell2DEdges x NumberCell2D

};

}
