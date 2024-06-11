#include <iostream>
#include "GeometryDFN.hpp"
#include "Utils.hpp"
#include "Utils2.hpp"
#include "fstream"
#include "Paraview/src/UCDUtilities.hpp"
#include <filesystem>

using namespace std;
using namespace Eigen;
using namespace DFN_Library;
using namespace DFN_PolygonalLibrary;


int main() {

    Struttura_DFN DFN;

    // Popola la struttura DFN con i dati necessari
    string filename = "./FR10_data.txt";
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
    string OutputNameTraccia = "./results/Tracce_FR10.txt";
    string OutputNameFrattura = "./results/Fratture_FR10.txt";
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
    string OutputNameMesh = "./resultsParte2/Tracce_FR10";
    if (!OutputPolygonalMesh(DFN, OutputNameMesh)){
        return 6;
    }

    // Paraview
    unsigned int id_frattura = 0;
    vector<list<Vector3d>> insiemePoligoni = trovaPoligoniTotali(id_frattura, DFN);
    PolygonalMesh Mesh = calcolaCelle0D(insiemePoligoni);
    calcolaCelle1D2D(insiemePoligoni, Mesh);
    vector<vector<unsigned int>> sottoPoligoni = Mesh.Cell2DVertices;
    MatrixXd MatricePunti;
    for (unsigned int i = 0; i < Mesh.NumberCell0D; i++)
    {
        MatricePunti.col(i) = Mesh.Cell0DCoordinates[i];
    }
    Gedim::UCDUtilities exporter;
    ofstream ofname;
    ofname.open("./Poligono0_FR10.inp");
    exporter.ExportPolygons("./Poligono0_FR10.inp", MatricePunti, sottoPoligoni, {}, {}, {});
    ofname.close();

    return 0;
}
