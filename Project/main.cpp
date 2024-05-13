#include <iostream>
#include "GeometryDFN.hpp"
#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using namespace DFN_Library;


int main() {

    Struttura_DFN DFN;
    // Popola la struttura DFN con i dati necessari

    string filename = "./FR3_data.txt"; // Il nome del file in cui esportare i dati
    if(!ImportFratture(filename, DFN)){
        return 1;
    }

    // std::vector<Vector3d> poligono;

    // // Aggiungi i vertici del quadrato al vettore
    // Vector3d vertice1 = {0.0, 0.0, 0.0};
    // Vector3d vertice2 = {1.0, 0.0, 0.0};
    // Vector3d vertice3 = {1.0, 1.0, 0.0};
    // Vector3d vertice4 = {0.0, 1.0, 0.0};

    // // Aggiungi i vertici al vettore
    // poligono.push_back(vertice1);
    // poligono.push_back(vertice2);
    // poligono.push_back(vertice3);
    // poligono.push_back(vertice4);

    // // Vector3d centroide = calcolaCentroide(poligono);
    // // cout << "Centroide: (" << centroide[0] << ";" << centroide[1] << ";" << centroide[2] << endl;

    // std::vector<Vector3d> poligono2;

    // // Aggiungi i vertici del quadrato al vettore
    // Vector3d vertice5 = {3.0, 3.0, 0.0};
    // Vector3d vertice6 = {3.5, 3.0, 0.0};
    // Vector3d vertice7 = {3.5, 3.5, 0.0};
    // Vector3d vertice8 = {0.0, 3.5, 0.0};

    // // Aggiungi i vertici al vettore
    // poligono2.push_back(vertice5);
    // poligono2.push_back(vertice6);
    // poligono2.push_back(vertice7);
    // poligono2.push_back(vertice8);

    // // Vector3d centroide2 = calcolaCentroide(poligono2);

    // bool risposta = possibiliTracce(poligono, poligono2, 1e-16);
    // cout << risposta << endl;
    // return 0;
}
