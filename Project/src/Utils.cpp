#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Eigen/Eigen"
#include "cmath"

using namespace std;
using namespace Eigen;

namespace DFN_Library{


bool funzioneMadre(const string &filename, Struttura_DFN& DFN)
{
    if (ImportFratture(filename, DFN))
    {
        unsigned int contatoreTraccia = 0;

        for (unsigned int i = 0; i < DFN.numFratture; i++) // Codice id della prima frattura
        {
            for (unsigned int j = i+1; j < DFN.numFratture; j++) // Codice id della seconda frattura
            {
                if (possibiliTracce(DFN.coordinateVertici[i], DFN.coordinateVertici[j], 1e-5))
                {
                    Vector3d normaleFrattura_i = normalePoligono(DFN.coordinateVertici[i]);
                    Vector3d normaleFrattura_j = normalePoligono(DFN.coordinateVertici[j]);
                    Vector3d versoreTangente = normaleFrattura_i.cross(normaleFrattura_j);
                    Matrix3d MatriceA;
                    MatriceA.row(0) = normaleFrattura_i.transpose();
                    MatriceA.row(1) = normaleFrattura_j.transpose();
                    MatriceA.row(2) = versoreTangente.transpose();

                    double bi = normaleFrattura_i.transpose() * calcolaCentroide(DFN.coordinateVertici[i]);
                    double bj = normaleFrattura_j.transpose() * calcolaCentroide(DFN.coordinateVertici[j]);
                    Vector3d b = {bi, bj, 0};

                    if (MatriceA.determinant() < 1e-9)
                    {
                        // non c'è una soluzione/intersezione tra i piani
                        continue;
                    }
                    Vector3d puntoTraccia = MatriceA.fullPivLu().solve(b);

                    unsigned int numPuntiTracce = 0;
                    for(unsigned int k = 0; k < DFN.numVertici[i]; k++)
                    {
                        Vector3d vertice1 = DFN.coordinateVertici[i][k];
                        Vector3d vertice2 = DFN.coordinateVertici[i][0];
                        if (k != DFN.numVertici[i] - 1)
                        {
                            Vector3d vertice2 = DFN.coordinateVertici[i][k+1];
                        }

                        Matrix3d MatriceA1; // deve dicentare matrice 3x2
                        MatriceA.row(0) = (vertice2-vertice1);
                    }


                }
            }

        }
    }

    return false;
}

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

Vector3d calcolaCentroide(vector<Vector3d>& poligono)
{
    //const auto& vertici = DFN.coordinateVertici[idFrattura];   //richiamo delle coordinate dei vertici relativia frattura di Id dato
    Vector3d centroide = Vector3d::Zero();    //inizializzo a zero il centroide
    int n = poligono.size();
    for (const auto& vertice : poligono) {
        centroide += vertice;
    }
    centroide /= n;     //calcolo del centroide come somma dei vertici diviso il loro totale
    return centroide;
}

bool possibiliTracce(vector<Vector3d>& poligono1, vector<Vector3d>& poligono2, double tolleranza)
{
    // Controllo sulla distanza delle fratture:
    Vector3d centroide1 = calcolaCentroide(poligono1);
    Vector3d centroide2 = calcolaCentroide(poligono2);

    double distanzaMax1 = 0;
    for (unsigned int j = 0; j < poligono1.size(); j++) // Modificato per scorrere tutti i vertici di poligono1
    {
        double distanza = (poligono1[j] - centroide1).squaredNorm();
        if (distanza > distanzaMax1)
        {
            distanzaMax1 = distanza; // Corretto per assegnare il valore massimo trovato
        }
    }

    double distanzaMax2 = 0;
    for (unsigned int j = 0; j < poligono2.size(); j++) // Modificato per scorrere tutti i vertici di poligono1
    {
        double distanza = (poligono2[j] - centroide2).squaredNorm();
        if (distanza > distanzaMax2)
        {
            distanzaMax2 = distanza; // Corretto per assegnare il valore massimo trovato
        }
    }

    double distanzaCentroidi = (centroide1 - centroide2).squaredNorm();

    // // Controllo sulla distanza delle fratture:
    // double distanzaCentroidi = 0;
    // for (unsigned int j = 0; j < 3; j++){
    //     distanzaCentroidi += pow((centroide1[j] - centroide2[j]), 2);
    // }
    // double distanzaMax1 = 0;
    // for (unsigned int j = 0; j < poligono1.size(); j++) // Modificato per scorrere tutti i vertici di poligono1
    // {
    //     double distanza = 0;
    //     for (unsigned int k = 0; k < 3; k++){
    //         distanza += pow((poligono1[j][k] - centroide1[k]), 2);
    //     }
    //     if (distanza > distanzaMax1)
    //     {
    //         distanzaMax1 = distanza; // Corretto per assegnare il valore massimo trovato
    //     }
    // }

    // double distanzaMax2 = 0;
    // for (unsigned int j = 0; j < poligono2.size(); j++) // Modificato per scorrere tutti i vertici di poligono1
    // {
    //     double distanza = 0;
    //     for (unsigned int k = 0; k < 3; k++){
    //         distanza += pow((poligono2[j][k] - centroide2[k]), 2);
    //     }
    //     if (distanza > distanzaMax2)
    //     {
    //         distanzaMax2 = distanza; // Corretto per assegnare il valore massimo trovato
    //     }
    // }

    // // Controllo sulla distanza delle fratture:
    // double distanzaCentroidi = 0;
    // for (unsigned int j = 0; j < 3; j++){
    //     distanzaCentroidi += pow((centroide1[j] - centroide2[j]), 2);
    // }


    if ((distanzaMax1 + distanzaMax2 + tolleranza) > distanzaCentroidi)
    {
        // Le fratture si possono incontrare
        return true; // Modificato per uscire dalla funzione in caso di fratture vicine a sufficienza
        // Se viene superato questo controllo posso procedere con il sistema lineare che mi calcoli le fratture
    }
    return false; // Modificato per indicare se le fratture sono troppo distanti
}

Vector3d normalePoligono(vector<Vector3d>& poligono)
{
    Vector3d punto0 = poligono[0];
    Vector3d punto1 = poligono[1];
    Vector3d punto2 = poligono[2];

    Vector3d vettore1_0 = punto1 - punto0; // vettore che parte da p1 e arriva in p0
    Vector3d vettore1_2 = punto1 - punto2; // vettore che parte da p1 e arriva in p2

    Vector3d normale = (vettore1_0.cross(vettore1_2))/(vettore1_0.norm() * vettore1_2.norm());
    return normale;
}


}

