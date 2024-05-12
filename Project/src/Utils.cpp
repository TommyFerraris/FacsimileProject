#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace DFN_Library{


bool ImportFratture(const string &filename,
                    Struttura_DFN& DFN)
{

    ifstream file;
    file.open(filename);

    if(file.fail())
        return false;

    string line;
    //leggo la riga senza considerarla perchè commento
    getline(file, line);
    getline(file,line);


    istringstream convertNumber;
    convertNumber.str(line);
    convertNumber >> DFN.numFratture;

    DFN.Id_Fratture.reserve(DFN.numFratture);
    DFN.numVertici.reserve(DFN.numFratture);

    for(unsigned int i = 0; i < DFN.numFratture; i++)
    {
        getline(file, line); // ignoro commento
        getline(file, line);
        replace(line.begin(), line.end(), ';', ' ');
        istringstream converter(line);

        unsigned int id;
        unsigned int numeroVertici;
        converter >> id >> numeroVertici;

        DFN.Id_Fratture.push_back(id);
        DFN.numVertici.push_back(numeroVertici);

        getline(file,line); //ignoro commento
        vector<Vector3d> matriceVertici;
        matriceVertici.resize(numeroVertici);


        getline(file, line);
        vector<double> vettorecoordinatex;
        stringstream ss(line);
        string token;
        while (getline(ss, token, ';')) {
            vettorecoordinatex.push_back(stod(token));
        }


        for (unsigned int k=0; k < numeroVertici; k++)
        {
            matriceVertici[k][0] = vettorecoordinatex[k];

        }

        getline(file, line);
        vector<double> vettorecoordinatey;
        stringstream ss2(line);
        string token2;
        while (getline(ss2, token2, ';')) {
            vettorecoordinatey.push_back(stod(token2));
        }


        for (unsigned int k=0; k < numeroVertici; k++)
        {
            matriceVertici[k][1] = vettorecoordinatey[k];

        }

        getline(file, line);
        vector<double> vettorecoordinatez;
        stringstream ss3(line);
        string token3;
        while (getline(ss3, token3, ';')) {
            vettorecoordinatez.push_back(stod(token3));
        }


        for (unsigned int k=0; k < numeroVertici; k++)
        {
            matriceVertici[k][2] = vettorecoordinatez[k];

        }

        DFN.coordinateVertici[id] = matriceVertici;
    }
    file.close();

    for(unsigned int i = 0; i < DFN.numFratture; i++)
    {
        // Codice precedente per leggere le informazioni della frattura

        // Stampa dei valori letti per la frattura corrente
        cout << "Frattura ID: " << DFN.Id_Fratture[i] << endl;
        cout << "Numero di vertici: " << DFN.numVertici[i] << endl;
        cout << "Coordinate dei vertici:" << endl;
        for (unsigned int k = 0; k < DFN.numVertici[i]; k++)
        {
            Vector3d vertice = DFN.coordinateVertici[DFN.Id_Fratture[i]][k];
            cout << "Vertice " << k << ": " << vertice.transpose() << endl;
        }
        cout << endl;
    }

    return true;
}

Vector3d calcolaCentroide(vector<Vector3d>& vertici)
{
    //const auto& vertici = DFN.coordinateVertici[idFrattura];   //richiamo delle coordinate dei vertici relativia frattura di Id dato
    Vector3d centroide = Vector3d::Zero();    //inizializzo a zero il centroide
    int n = vertici.size();
    for (const auto& vertex : vertici) {
        centroide += vertex;
    }
    centroide /= n;                           //calcolo del centroide come somma dei vertici diviso il loro totale
    return centroide;
}

bool CalcolaTracce(Struttura_DFN& DFN)
{
    for (unsigned int i = 0; i < DFN.numFratture; i++) // Modificato per evitare l'accesso out-of-bounds
    {
        const auto& verticiFrattura1 = DFN.coordinateVertici[DFN.Id_Fratture[i]];
        const auto& verticiFrattura2 = DFN.coordinateVertici[DFN.Id_Fratture[i+1]];

        // Controllo sulla distanza delle fratture:
        Vector3d centroide1 = calcolaCentroide(DFN, DFN.Id_Fratture[i]);
        Vector3d centroide2 = calcolaCentroide(DFN, DFN.Id_Fratture[i+1]);

        double distanzaMax1 = 0;
        for (unsigned int j = 0; j < verticiFrattura1.size(); j++) // Modificato per scorrere tutti i vertici di verticiFrattura1
        {
            double distanza = (verticiFrattura1[j] - centroide1).norm();
            if (distanza > distanzaMax1)
            {
                distanzaMax1 = distanza; // Corretto per assegnare il valore massimo trovato
            }
        }

        double distanzaMax2 = 0;
        for (unsigned int j = 0; j < verticiFrattura2.size(); j++) // Modificato per scorrere tutti i vertici di verticiFrattura2
        {
            double distanza = (verticiFrattura2[j] - centroide2).norm();
            if (distanza > distanzaMax2)
            {
                distanzaMax2 = distanza; // Corretto per assegnare il valore massimo trovato
            }
        }

        // Controllo sulla distanza delle fratture:
        double distanzaCentroidi = (centroide1 - centroide2).norm();

        if ((distanzaMax1 + distanzaMax2) > distanzaCentroidi)
        {
            // Le fratture sono troppo lontane
            return false; // Modificato per uscire dalla funzione in caso di fratture troppo distanti
        }
        // Se viene superato questo controllo posso procedere con il sistema lineare che mi calcoli le fratture:
    }
    return true; // Modificato per indicare che tutte le fratture sono state processate con successo
}


}

