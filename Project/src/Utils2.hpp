#pragma once

#include "PolygonalMesh.hpp"
#include "GeometryDFN.hpp"

using namespace std;
using namespace DFN_Library;


namespace DFN_PolygonalLibrary {
// Funzione che dato un poligono e i due punti di intersezione con la traccia, me lo divide
array<vector<Vector3d>,2> dividiPoligono(vector<Vector3d>& poligono, Vector3d& p1, Vector3d& p2);
// Funzione che dato un poligono e le sue tracce, me lo suddivide nei poligoni che appartengono alle celle 2D della mesh
pair<vector<vector<Vector3d>>, PolygonalMesh> trovaPoligoniTotali(unsigned int& Idpoligono, Struttura_DFN& DFN);
// Funzione che ci salvi in una mesh le celle1D
void calcolaCelle1D (vector<vector<Vector3d>>& insiemePoligoni, PolygonalMesh& Mesh);
// Funzione che mi calcoli il fattoriale, costruita in modo ricorsivo
unsigned long long fattoriale (unsigned int n);
}
