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

    return 0;
}

