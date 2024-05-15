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

    Vector3d centroide0 = calcolaCentroide(DFN.coordinateVertici[0]);
    Vector3d centroide1 = calcolaCentroide(DFN.coordinateVertici[1]);
    Vector3d centroide2 = calcolaCentroide(DFN.coordinateVertici[2]);
    cout << "Cenrtroide 0: (" << centroide0[0] << ";" << centroide0[1] << ";" << centroide0[2] << endl;
    cout << "Cenrtroide 1: (" << centroide1[0] << ";" << centroide1[1] << ";" << centroide1[2] << endl;
    cout << "Cenrtroide 2: (" << centroide2[0] << ";" << centroide2[1] << ";" << centroide2[2] << endl;

    bool risposta01 = possibiliTracce(DFN.coordinateVertici[0], DFN.coordinateVertici[1], 5e-1);
    bool risposta02 = possibiliTracce(DFN.coordinateVertici[0], DFN.coordinateVertici[2], 5e-1);
    bool risposta12 = possibiliTracce(DFN.coordinateVertici[1], DFN.coordinateVertici[2], 5e-1);
    cout << risposta01 << "; " << risposta02 << "; " << risposta12 << endl;

    Vector3d normale = normalePoligono(DFN.coordinateVertici[0]);
    cout << "Normale: (" << normale[0] << ";" << normale[1] << ";" << normale[2] << endl;

    // return 0;
}
