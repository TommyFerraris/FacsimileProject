#pragma once

#include "PolygonalMesh.hpp"
#include "GeometryDFN.hpp"

using namespace std;


namespace DFN_PolygonalLibrary {
array<vector<Vector3d>,2> dividiPoligono(vector<Vector3d>& poligono, Vector3d& p1, Vector3d& p2);
}
