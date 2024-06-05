#include <iostream>
#include "GeometryDFN.hpp"
#include "Utils.hpp"
#include "Utils2.hpp"
#include <filesystem>

using namespace std;
using namespace Eigen;
using namespace DFN_Library;
using namespace DFN_PolygonalLibrary;


int main() {

    Struttura_DFN DFN;

    // Popola la struttura DFN con i dati necessari
    string filename = "./FR3_data.txt";
    if(!ImportFratture(filename, DFN)){
        return 1;
    }

    // Funzioni che lavorano sulla struttura DFN
    calcolaTracce(DFN);
    calcolaTipologiaTracce(DFN);
    calcolaLunghezzaTracce(DFN);

    // Funzioni che mi salvino i risultati in file .txt
    string directory = "results";
    if (!filesystem::exists(directory)) {
        if (!filesystem::create_directory(directory)) {
            std::cerr << "Errore nella creazione della cartella: " << directory << std::endl;
            return 2;
        }
    }
    string OutputNameTraccia = "./results/Tracce_FR3.txt";
    string OutputNameFrattura = "./results/Fratture_FR3.txt";
    if(!OutputTracce(DFN, OutputNameTraccia)){
        return 3;
    }
    if(!OutputFratture(DFN, OutputNameFrattura)){
        return 4;
    }

    // Funzioni che mi salvino i risultati in file .txt
    string directory2 = "resultsParte2";
    if (!filesystem::exists(directory2)) {
        if (!filesystem::create_directory(directory2)) {
            std::cerr << "Errore nella creazione della cartella: " << directory2 << std::endl;
            return 5;
        }
    }
    string OutputNameMesh = "./resultsParte2/Tracce_FR3";
    if (!OutputPolygonalMesh(DFN, OutputNameMesh)){
        return 6;
    }

    return 0;
}
