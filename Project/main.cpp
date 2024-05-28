#include <iostream>
#include "GeometryDFN.hpp"
#include "Utils.hpp"
#include <filesystem>

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

    // Funzioni che lavorano sulla struttura DFN
    calcolaTracce(DFN);
    calcolaTipologiaTracce(DFN);
    calcolaLunghezzaTracce(DFN);

    // Funzioni che mi salvino i risultati in file .csv
    string directory = "results";
    if (!filesystem::exists(directory)) {
        if (!filesystem::create_directory(directory)) {
            std::cerr << "Errore nella creazione della cartella: " << directory << std::endl;
            return 4;
        }
    }
    string OutputNameTraccia = "./results/Tracce_FR3.csv";
    string OutputNameFrattura = "./results/Fratture_FR3.csv";
    if(!OutputTracce(DFN, OutputNameTraccia)){
        return 2;
    }
    if(!OutputFratture(DFN, OutputNameFrattura)){
        return 3;
    }

    // for (unsigned int i = 0; i < DFN.numTracce; i++)
    // {
    //     cout << "Id traccia: " << i << ", lunghezza: " << DFN.lunghezzaTraccia[i] << endl;
    // }


    // for (unsigned int i = 0; i < DFN.numFratture; i++)
    // {
    //     if (DFN.tipoTraccia.count(i)>0)
    //     {
    //         cout << "Per la frattura " << i << " ci sono tali tracce: " << endl;
    //         for(unsigned int j = 0; j < DFN.tipoTraccia[i].size(); j++)
    //         {
    //             int id_traccia = DFN.tipoTraccia[i][j][0];
    //             int tipotraccia = DFN.tipoTraccia[i][j][1];
    //             cout << "per la traccia " << id_traccia << " la sua tipologia equivale a: " << tipotraccia << endl;
    //         }
    //     }
    // }

    // for(unsigned int i = 0; i < DFN.numTracce; i++)
    // {
    //     cout << "Id traccia: " << DFN.Id_Traccia[i] << endl;
    //     cout << "Id fratture intersecanti: " << DFN.Id_Fratture_Intersecanti[i][0] << "; " << DFN.Id_Fratture_Intersecanti[i][1];
    //     cout << "Coordinate traccia: " << DFN.coordinateTraccia[i][0][0] << ", " << DFN.coordinateTraccia[i][0][1] << ", " <<
    //         DFN.coordinateTraccia[i][0][2] << ", secondo punto: " << DFN.coordinateTraccia[i][1][0] << ", " <<
    //         DFN.coordinateTraccia[i][1][1] << ", " << DFN.coordinateTraccia[i][1][2] << ", " << endl;
    // }

    // Vector3d centroide0 = calcolaCentroide(DFN.coordinateVertici[0]);
    // Vector3d centroide1 = calcolaCentroide(DFN.coordinateVertici[1]);
    // Vector3d centroide2 = calcolaCentroide(DFN.coordinateVertici[2]);
    // cout << "Cenrtroide 0: (" << centroide0[0] << ";" << centroide0[1] << ";" << centroide0[2] << endl;
    // cout << "Cenrtroide 1: (" << centroide1[0] << ";" << centroide1[1] << ";" << centroide1[2] << endl;
    // cout << "Cenrtroide 2: (" << centroide2[0] << ";" << centroide2[1] << ";" << centroide2[2] << endl;

    // bool risposta01 = possibiliTracce(DFN.coordinateVertici[0], DFN.coordinateVertici[1], 5e-1);
    // bool risposta02 = possibiliTracce(DFN.coordinateVertici[0], DFN.coordinateVertici[2], 5e-1);
    // bool risposta12 = possibiliTracce(DFN.coordinateVertici[1], DFN.coordinateVertici[2], 5e-1);
    // cout << risposta01 << "; " << risposta02 << "; " << risposta12 << endl;

    // Vector3d normale = normalePoligono(DFN.coordinateVertici[0]);
    // cout << "Normale: (" << normale[0] << ";" << normale[1] << ";" << normale[2] << endl;

    return 0;
}
