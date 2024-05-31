#pragma once

#include "PolygonalMesh.hpp"
#include "GeometryDFN.hpp"

using namespace std;
using namespace DFN_Library;


namespace DFN_PolygonalLibrary {
// Funzione che dato un poligono e i due punti di intersezione con la traccia, me lo divide
array<vector<Vector3d>,2> dividiPoligono(vector<Vector3d>& poligono, Vector3d& p1, Vector3d& p2);
// Funzione che cpntrolla se la traccia divide il poligono o meno
bool trovaIntersezioneTracciaPoligono(vector<Vector3d>& poligono, Vector3d& p1, Vector3d& p2);
vector<vector<Vector3d>> trovaPoligoniTotali(unsigned int& Idpoligono, Struttura_DFN& DFN);
}
