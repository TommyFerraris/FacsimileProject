#include "Utils2.hpp"
#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Eigen"
#include "cmath"
#include "algorithm"

using namespace std;
using namespace Eigen;
using namespace DFN_Library;

namespace DFN_PolygonalLibrary{

array<list<Vector3d>, 2> dividiPoligono(list<Vector3d>& poligono, const Vector3d& p1, const Vector3d& p2) {
    auto p1It = poligono.end();
    auto p2It = poligono.end();

    // Find positions of p1 and p2
    for (auto it = poligono.begin(); it != poligono.end(); ++it) {
        auto nextIt = next(it);
        if (nextIt == poligono.end()) nextIt = poligono.begin();

        if (puntoInSegmento(*it, *nextIt, p1)) {
            p1It = it;
        }
        if (puntoInSegmento(*it, *nextIt, p2)) {
            p2It = it;
        }
    }
    // Create two new polygons
    list<Vector3d> poligono1, poligono2;

    // Fill poligono1
    poligono1.push_back(p2);
    for (auto it = (p2It); it != (p1It);) {
        it = next(it);
        if (it == poligono.end()) it = poligono.begin();
        poligono1.push_back(*it);
    }
    poligono1.push_back(p1);

    // Fill poligono2
    poligono2.push_back(p1);
    for (auto it = p1It; it != p2It;) {
        it = next(it);
        if (it == poligono.end()) it = poligono.begin();
        poligono2.push_back(*it);
    }
    poligono2.push_back(p2);

    return {poligono1, poligono2};
}

bool puntoInternoPoligonoLista(const Vector3d& punto, const std::list<Vector3d>& poligono)
{
    auto it = poligono.begin();
    Vector3d primoVertice = *it;  // Primo vertice del poligono
    ++it;  // Avanza al secondo vertice

    for (auto prev = it, next = std::next(it); next != poligono.end(); ++prev, ++next)
    {
        if (puntointriangolo(punto, primoVertice, *prev, *next))
        {
            return true;
        }
    }

    return false;
}
vector<list<Vector3d>> trovaPoligoniTotali(unsigned int& Idpoligono, Struttura_DFN& DFN) {
    vector<list<Vector3d>> vettorePoligoni;

    // Convert vector<Vector3d> to list<Vector3d>
    list<Vector3d> poligono(DFN.coordinateVertici[Idpoligono].begin(), DFN.coordinateVertici[Idpoligono].end());
    vettorePoligoni.push_back(poligono);

    for (unsigned int i = 0; i < DFN.tipoTraccia[Idpoligono].size(); i++) {
        unsigned int grandezza = vettorePoligoni.size();
        Vector3d OrigineTraccia = DFN.coordinateTraccia[DFN.tipoTraccia[Idpoligono][i][0]][0];
        Vector3d FineTraccia = DFN.coordinateTraccia[DFN.tipoTraccia[Idpoligono][i][0]][1];

        for (unsigned int j = 0; j < grandezza; j++) {
            list<Vector3d>& currentPolygon = vettorePoligoni[j];
            vector<Vector3d> puntiIntersezione;
            unsigned int contatore = 0;
            if (puntoInternoPoligonoLista(OrigineTraccia, currentPolygon) && puntoInternoPoligonoLista(FineTraccia, currentPolygon))
            {
                contatore += 1;
            }
            for (auto it = currentPolygon.begin(); it != currentPolygon.end(); ++it) {
                auto nextIt = next(it);
                if (nextIt == currentPolygon.end()) nextIt = currentPolygon.begin();

                Vector3d vertice1 = *it;
                Vector3d vertice2 = *nextIt;

                Matrix<double, 3, 2> MatriceA1;
                MatriceA1.col(0) = (vertice2 - vertice1);
                MatriceA1.col(1) = (FineTraccia - OrigineTraccia);

                if ((vertice2 - vertice1).cross(FineTraccia - OrigineTraccia).squaredNorm() < 1e-12) {
                    continue; // Lines are parallel
                }

                Vector3d b1 = (OrigineTraccia - vertice1);
                Vector2d alphaBeta = MatriceA1.colPivHouseholderQr().solve(b1);

                if (alphaBeta[0] >= -1e-09 && alphaBeta[0] <= 1 + 1e-09 &&
                    abs(alphaBeta[1]) >= -1e-09 && abs(alphaBeta[1]) <= 1 + 1e-09) {
                    Vector3d Intersezione = vertice1 + alphaBeta[0] * (vertice2 - vertice1);
                    puntiIntersezione.push_back(Intersezione);
                    contatore += 1;
                }
                else if (alphaBeta[0] >= -1e-09 && alphaBeta[0] <= 1 + 1e-09)
                {
                    Vector3d Intersezione = vertice1 + alphaBeta[0] * (vertice2 - vertice1);
                    puntiIntersezione.push_back(Intersezione);
                }
            }

            if (contatore != 0 && puntiIntersezione.size() == 2) {
                array<list<Vector3d>, 2> poligoniNuovi = dividiPoligono(currentPolygon, puntiIntersezione[0], puntiIntersezione[1]);
                replace(vettorePoligoni.begin(), vettorePoligoni.end(), currentPolygon, poligoniNuovi[0]);
                vettorePoligoni.push_back(poligoniNuovi[1]);
            }
        }
    }

    return vettorePoligoni;
}


PolygonalMesh calcolaCelle0D(vector<list<Vector3d>>& insiemePoligoni)
{
    PolygonalMesh Mesh;
    unsigned int riservaPosti = 0;
    for(unsigned int i = 0; i < insiemePoligoni.size(); i++) {riservaPosti += insiemePoligoni[i].size();}
    Mesh.Cell0DId.reserve(riservaPosti);
    Mesh.Cell0DCoordinates.reserve(riservaPosti);
    unsigned int Id = 0;

    for (unsigned int i = 0; i < insiemePoligoni.size(); i++) {
        for (auto it = insiemePoligoni[i].begin(); it != insiemePoligoni[i].end(); ++it) {
            Vector3d vertice = *it;
            if (find(Mesh.Cell0DCoordinates.begin(), Mesh.Cell0DCoordinates.end(), vertice) == Mesh.Cell0DCoordinates.end())
            {
                Mesh.NumberCell0D += 1;
                Mesh.Cell0DId.push_back(Id);
                Mesh.Cell0DCoordinates.push_back(vertice);
                Id += 1;
            }
        }
    }
    Mesh.Cell0DCoordinates.shrink_to_fit();
    Mesh.Cell0DId.shrink_to_fit();
    return Mesh;
}


void calcolaCelle1D2D(vector<list<Vector3d>>& insiemePoligoni, PolygonalMesh& Mesh) {
    unsigned int Id = 0;
    Mesh.Cell2DId.reserve(insiemePoligoni.size());
    Mesh.Cell2DNumVertices.reserve(insiemePoligoni.size());
    Mesh.Cell2DVertices.reserve(insiemePoligoni.size());
    Mesh.Cell2DNumEdges.reserve(insiemePoligoni.size());
    Mesh.Cell2DEdges.reserve(insiemePoligoni.size());
    unsigned int Id2D = 0;
    Mesh.NumberCell2D = insiemePoligoni.size();

    for (unsigned int i = 0; i < insiemePoligoni.size(); i++) {
        Mesh.Cell2DId.push_back(Id2D);
        Mesh.Cell2DNumVertices.push_back(insiemePoligoni[i].size());
        Mesh.Cell2DNumEdges.push_back(insiemePoligoni[i].size());
        Id2D += 1;
        vector<unsigned int> Vertici;
        vector<unsigned int> Lati;
        for (auto it = insiemePoligoni[i].begin(); it != insiemePoligoni[i].end(); ++it) {
            Vector3d vertice1 = *it;
            auto IdV1 = find(Mesh.Cell0DCoordinates.begin(), Mesh.Cell0DCoordinates.end(), vertice1);
            int Id1 = distance(Mesh.Cell0DCoordinates.begin(), IdV1);
            Vertici.push_back(Id1);

            auto nextIt = next(it);
            if (nextIt == insiemePoligoni[i].end()) {
                nextIt = insiemePoligoni[i].begin();
            }
            Vector3d vertice2 = *nextIt;
            auto IdV2 = find(Mesh.Cell0DCoordinates.begin(), Mesh.Cell0DCoordinates.end(), vertice2);
            int Id2 = distance(Mesh.Cell0DCoordinates.begin(), IdV2);

            Vector2i prova1 = {Id1, Id2};
            Vector2i prova2 = {Id2, Id1};

            if (find(Mesh.Cell1DVertices.begin(), Mesh.Cell1DVertices.end(), prova1) == Mesh.Cell1DVertices.end() &&
                find(Mesh.Cell1DVertices.begin(), Mesh.Cell1DVertices.end(), prova2) == Mesh.Cell1DVertices.end()) {
                Mesh.NumberCell1D += 1;
                Mesh.Cell1DId.push_back(Id);
                Mesh.Cell1DVertices.push_back(prova1);
                Lati.push_back(Id);
                Id += 1;
            }
            else
            {
                auto IdE1 = find(Mesh.Cell1DVertices.begin(), Mesh.Cell1DVertices.end(), prova1);
                auto IdE2 = find(Mesh.Cell1DVertices.begin(), Mesh.Cell1DVertices.end(), prova2);
                if (IdE1 != Mesh.Cell1DVertices.end())
                {
                    int IdE = distance(Mesh.Cell1DVertices.begin(), IdE1);
                    Lati.push_back(IdE);
                }
                else
                {
                    int IdE = distance(Mesh.Cell1DVertices.begin(), IdE2);
                    Lati.push_back(IdE);
                }
            }
        }
        Mesh.Cell2DVertices.push_back(Vertici);
        Mesh.Cell2DEdges.push_back(Lati);
    }
}


bool OutputPolygonalMesh(Struttura_DFN& DFN, const string& fileOutput)
{
    ofstream file1;
    file1.open(fileOutput+"Celle0D.txt");
    if (file1.fail())
    {
        cerr << "Non riesco ad aprire il file " << fileOutput << "Celle0D.txt" << endl;
        return false;
    }
    ofstream file2;
    file2.open(fileOutput+"Celle1D.txt");
    if (file2.fail())
    {
        cerr << "Non riesco ad aprire il file " << fileOutput << "Celle1D.txt" << endl;
        file1.close();
        return false;
    }
    ofstream file3;
    file3.open(fileOutput+"Celle2D.txt");
    if (file3.fail())
    {
        cerr << "Non riesco ad aprire il file " << fileOutput << "Celle2D.txt" << endl;
        file1.close();
        file2.close();
        return false;
    }

    for (unsigned int i = 0; i < DFN.Id_FrattureConTraccia.size(); i++)
    {
        unsigned int Id_Frattura = DFN.Id_FrattureConTraccia[i];
        vector<list<Vector3d>> insiemePoligoni = trovaPoligoniTotali(Id_Frattura, DFN);
        PolygonalMesh Mesh = calcolaCelle0D(insiemePoligoni);
        calcolaCelle1D2D(insiemePoligoni, Mesh);

        file1 << "# Id Frattura" << endl;
        file1 << Id_Frattura << endl;
        file1 << "# Id; X; Y; Z" << endl;
        file1.precision(16);
        file1 << scientific;
        for (unsigned int k = 0; k < Mesh.NumberCell0D; k++)
        {
            file1 << Mesh.Cell0DId[k] << "; " << Mesh.Cell0DCoordinates[k][0] << "; " << Mesh.Cell0DCoordinates[k][1] << "; " << Mesh.Cell0DCoordinates[k][2] << endl;
        }

        file2 << "# Id Frattura" << endl;
        file2 << Id_Frattura << endl;
        file2 << "#Id; Origin; End" << endl;
        for (unsigned int k = 0; k < Mesh.NumberCell1D; k++)
        {
            file2 << Mesh.Cell1DId[k] << "; " << Mesh.Cell1DVertices[k][0] << "; " << Mesh.Cell1DVertices[k][1] << endl;
        }

        file3 << "# Id Frattura" << endl;
        file3 << Id_Frattura << endl;
        file3 << "# Id; NumVertices; Vertices; NumEdges; Edges" << endl;
        for (unsigned int k = 0; k < Mesh.NumberCell2D; k++)
        {
            file3 << Mesh.Cell2DId[k] << "; " << Mesh.Cell2DNumVertices[k];
            for (unsigned int j = 0; j < Mesh.Cell2DNumVertices[k]; j++)
            {
                file3 << "; " << Mesh.Cell2DVertices[k][j];
            }
            file3 << "; " << Mesh.Cell2DNumEdges[k];
            for (unsigned int j = 0; j < Mesh.Cell2DNumEdges[k]; j++)
            {
                file3 << "; " << Mesh.Cell2DEdges[k][j];
            }
            file3 << endl;
        }
        file1 << endl;
        file2 << endl;
        file3 << endl;
    }
    file1.close();
    file2.close();
    file3.close();
    return true;
}

}
