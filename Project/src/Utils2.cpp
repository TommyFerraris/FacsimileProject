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
        poligono1.push_back(poligono[p2Posizione%poligono.size()]);
    }
    poligono1.push_back(p1);

    vector<Vector3d> poligono2;
    poligono2.push_back(p1);
    while((p1P + 1)%poligono.size() != (p2P + 1)%poligono.size())
    {
        p1P += 1;
        poligono2.push_back(poligono[p1P%poligono.size()]);
    }
    poligono2.push_back(p2);

    array<vector<Vector3d>,2> poligoni = {poligono1, poligono2};
    return poligoni;
}



vector<vector<Vector3d>> trovaPoligoniTotali(unsigned int& Idpoligono, Struttura_DFN& DFN)
{
    vector<vector<Vector3d>> vettorePoligoni;
    vector<Vector3d>& poligono = DFN.coordinateVertici[Idpoligono];
    vettorePoligoni.push_back(poligono);
    for (unsigned int i = 0; i < DFN.tipoTraccia[Idpoligono].size(); i++)
    {
        unsigned int grandezza = vettorePoligoni.size();
        Vector3d OrigineTraccia = DFN.coordinateTraccia[DFN.tipoTraccia[Idpoligono][i][0]][0];
        Vector3d FineTraccia = DFN.coordinateTraccia[DFN.tipoTraccia[Idpoligono][i][0]][1];
        for (unsigned int j = 0; j < grandezza; j++)
        {
            vector<Vector3d> puntiIntersezione;
            unsigned int contatore = 0;
            for(unsigned int k = 0; k < poligono.size(); k++)
            {
                Vector3d vertice1 = vettorePoligoni[j][k];
                Vector3d vertice2 = vettorePoligoni[j][(k + 1) % poligono.size()];

                MatrixXd MatriceA1(3,2); // deve diventare matrice 3x2
                MatriceA1.col(0) = (vertice2-vertice1);
                MatriceA1.col(1) = (FineTraccia - OrigineTraccia);
                if ((vertice2-vertice1).cross(FineTraccia - OrigineTraccia).squaredNorm() < 1e-12)
                {
                    // le due rette sono parallele e quindi le escludo
                    continue;
                }
                Vector3d b1 = (OrigineTraccia - vertice1);

                Vector2d alphaBeta = MatriceA1.fullPivLu().solve(b1);
                if (alphaBeta[0] < 1 + 1e-09 && alphaBeta[0] > -1e-09 && alphaBeta[1] < 1 + 1e-09 && alphaBeta[1] > -1e-09)
                {
                    Vector3d Intersezione = vertice1 + (alphaBeta[0] * (vertice2 - vertice1));
                    puntiIntersezione.push_back(Intersezione);
                    contatore += 1;
                }
                else if (alphaBeta[0] < 1 + 1e-09 && alphaBeta[0] > -1e-09)
                {
                    Vector3d Intersezione = vertice1 + (alphaBeta[0] * (vertice2 - vertice1));
                    puntiIntersezione.push_back(Intersezione);
                }
            }
            if (contatore != 0 && puntiIntersezione.size() == 2)
            {
                array<vector<Vector3d>,2> poligoniNuovi= dividiPoligono(vettorePoligoni[j], puntiIntersezione[0], puntiIntersezione[1]);
                vettorePoligoni[j] = poligoniNuovi[0];
                vettorePoligoni.push_back(poligoniNuovi[1]);
            }
        }
    }
    return vettorePoligoni;
}
}
