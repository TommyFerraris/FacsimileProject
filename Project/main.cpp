#include <iostream>
#include "GeometryDFN.hpp"
#include "Utils.hpp"
#include "Utils2.hpp"
// #include "UCDUtilities.hpp"       //da scommentare se si vuole stampare in paraview
// #include "Triangolazione.hpp"     //da scommentare se si vuole stampare in paraview
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

    // Funzioni che mi salvano i risultati in file .txt
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

    // Funzioni che salvano i risultati in un file .txt
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

    // Paraview (per visualizzare esempi di alcune mesh)

    // unsigned int id_frattura = 1;
    // vector<list<Vector3d>> insiemePoligoni = trovaPoligoniTotali(id_frattura, DFN);
    // PolygonalMesh Mesh = calcolaCelle0D(insiemePoligoni);
    // calcolaCelle1D2D(insiemePoligoni, Mesh);
    // vector<vector<unsigned int>> sottoPoligoni = Mesh.Cell2DVertices;
    // MatrixXd MatricePunti (3, Mesh.NumberCell0D);
    // for (unsigned int i = 0; i < Mesh.NumberCell0D; i++)
    // {
    //     MatricePunti.col(i) = Mesh.Cell0DCoordinates[i];
    // }
    // GeometryLibrary::Polygons polygons;
    // polygons.VerticesCoordinates = MatricePunti;
    // polygons.listVertices = sottoPoligoni;
    // Gedim::UCDUtilities exporter;
    // std::vector<std::vector<unsigned int>> triangles;
    // Eigen::VectorXi materials;
    // polygons.GedimInterface(triangles, materials);
    // exporter.ExportPolygons("./Frattura50_1.inp",
    //                         polygons.VerticesCoordinates,
    //                         triangles,
    //                         {},
    //                         {},
    //                         materials);

    // VectorXi materials0D(Mesh.NumberCell0D);
    // for (unsigned int i = 0; i < Mesh.NumberCell0D; i++)
    // {
    //     materials0D(i) = Mesh.Cell0DId[i];
    // }
    // exporter.ExportPoints("./Frattura50_1_Celle0D.inp",
    //                       MatricePunti,
    //                       {},
    //                       materials0D);

    // MatrixXi lati(2, Mesh.NumberCell1D);
    // for (unsigned int i = 0; i < Mesh.NumberCell1D; i++)
    // {
    //     lati(0,i) = Mesh.Cell1DVertices[i][0];
    //     lati(1,i) = Mesh.Cell1DVertices[i][1];
    // }
    // VectorXi materials1D(Mesh.NumberCell1D);
    // for (unsigned int i = 0; i < Mesh.NumberCell1D; i++)
    // {
    //     materials1D(i) = Mesh.Cell1DId[i];
    // }
    // exporter.ExportSegments("./Frattura50_1_Celle1D.inp",
    //                         MatricePunti,
    //                         lati,
    //                         {},
    //                         {},
    //                         materials1D);

    return 0;
}
