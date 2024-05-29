#include "Utils2.hpp"
#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Eigen/Eigen"
#include "cmath"
#include "algorithm"

using namespace std;
using namespace Eigen;
using namespace DFN_Library;

namespace DFN_PolygonalLibrary{

array<vector<Vector3d>,2> dividiPoligono(vector<Vector3d>& poligono, Vector3d& p1, Vector3d& p2)
{
    unsigned int p1Posizione = 0;
    unsigned int p2Posizione = 0;
    for (unsigned int i = 0; i < poligono.size(); i++)
    {
        Vector3d V1 = poligono[i];
        Vector3d V2 = poligono[(i+1)%poligono.size()];
        if (puntoInSegmento(V1, V2, p1))
        {
            p1Posizione = i;
        }
        if (puntoInSegmento(V1, V2, p2))
        {
            p2Posizione = i;
        }
    }
    unsigned int p1P = p1Posizione;
    unsigned int p2P = p2Posizione;

    vector<Vector3d> poligono1;
    poligono1.push_back(p2);
    while((p2Posizione + 1)%poligono.size() != (p1Posizione + 1)%poligono.size())
    {
        p2Posizione += 1;
        poligono1.push_back(poligono[p2Posizione]);
    }
    poligono1.push_back(p1);

    vector<Vector3d> poligono2;
    poligono2.push_back(p1);
    while((p1P + 1)%poligono.size() != (p2P + 1)%poligono.size())
    {
        p1P += 1;
        poligono2.push_back(poligono[p1P]);
    }
    poligono2.push_back(p2);

    array<vector<Vector3d>,2> poligoni = {poligono1, poligono2};
    return poligoni;
}
}
