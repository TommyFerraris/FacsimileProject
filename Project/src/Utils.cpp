#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "Eigen/Eigen"
#include "cmath"
#include "algorithm"

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

    // for(unsigned int i = 0; i < DFN.numFratture; i++)
    // {
    //     // Codice precedente per leggere le informazioni della frattura

    //     // Stampa dei valori letti per la frattura corrente
    // //     cout << "Frattura ID: " << DFN.Id_Fratture[i] << endl;
    // //     cout << "Numero di vertici: " << DFN.numVertici[i] << endl;
    // //     cout << "Coordinate dei vertici:" << endl;
    // //     for (unsigned int k = 0; k < DFN.numVertici[i]; k++)
    // //     {
    // //         Vector3d vertice = DFN.coordinateVertici[DFN.Id_Fratture[i]][k];
    // //         cout << "Vertice " << k << ": " << vertice.transpose() << endl;
    // //     }
    // //     cout << endl;
    // }

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


void calcolaTracce(Struttura_DFN& DFN)
{
    unsigned int contatoreTraccia = 0;

    for (unsigned int i = 0; i < DFN.numFratture - 1; i++) // Codice id della prima frattura
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
                array<Vector3d,2> puntiTraccia;
                for(unsigned int k = 0; k < DFN.numVertici[i]; k++)
                {
                    Vector3d vertice1 = DFN.coordinateVertici[i][k];
                    Vector3d vertice2 = DFN.coordinateVertici[i][(k + 1) % DFN.numVertici[i]];

                    MatrixXd MatriceA1(3,2); // deve diventare matrice 3x2
                    MatriceA1.col(0) = (vertice2-vertice1);
                    MatriceA1.col(1) = -versoreTangente;
                    if ((vertice2-vertice1).cross(versoreTangente).squaredNorm() < 1e-09)
                    {
                        // le due rette sono parallele e quindi le escludo
                        continue;
                    }
                    Vector3d b1 = (puntoTraccia - vertice1);

                    Vector2d alphaBeta = MatriceA1.fullPivLu().solve(b1);
                    if (numPuntiTracce < 2)
                    {
                        Vector3d Intersezione = vertice1 + (alphaBeta[0] * (vertice2 - vertice1));
                        if (puntoInternoPoligono(Intersezione, DFN.coordinateVertici[j], normaleFrattura_j) && puntoInternoPoligono(Intersezione, DFN.coordinateVertici[i], normaleFrattura_i))
                        {
                            puntiTraccia[numPuntiTracce] = Intersezione;
                            numPuntiTracce += 1;
                        }
                    }
                }

                for(unsigned int k = 0; k < DFN.numVertici[j]; k++)
                {
                    Vector3d vertice1 = DFN.coordinateVertici[j][k];
                    Vector3d vertice2 = DFN.coordinateVertici[j][(k + 1) % DFN.numVertici[j]];
                    MatrixXd MatriceA1(3,2); // deve diventare matrice 3x2
                    MatriceA1.col(0) = (vertice2-vertice1);
                    MatriceA1.col(1) = -versoreTangente;
                    if ((vertice2-vertice1).cross(versoreTangente).squaredNorm() < 1e-09)
                    {
                        // le due rette sono parallele e quindi le escludo
                        continue;
                    }
                    Vector3d b1 = (puntoTraccia - vertice1);

                    Vector2d alphaBeta = MatriceA1.fullPivLu().solve(b1);
                    if (numPuntiTracce < 2)
                    {
                        Vector3d Intersezione = vertice1 + (alphaBeta[0] * (vertice2 - vertice1));
                        if (puntoInternoPoligono(Intersezione, DFN.coordinateVertici[i], normaleFrattura_i) && puntoInternoPoligono(Intersezione, DFN.coordinateVertici[j], normaleFrattura_j))
                        {
                            if (numPuntiTracce == 0)
                            {
                                puntiTraccia[numPuntiTracce] = Intersezione;
                                numPuntiTracce += 1;
                            }
                            else
                            {
                                puntiTraccia[numPuntiTracce] = Intersezione;
                                numPuntiTracce += 1;
                            }
                        }
                    }
                }

                DFN.numTracce += 1;
                DFN.Id_Traccia.push_back(contatoreTraccia);
                Vector2i fratture = {DFN.Id_Fratture[i], DFN.Id_Fratture[j]};
                DFN.Id_Fratture_Intersecanti.push_back(fratture);
                DFN.coordinateTraccia[contatoreTraccia] = puntiTraccia;
                contatoreTraccia += 1;

            }
        }
    }
}

//funzione calcolo del pun to interno con matrice di rotazione (NON VIENE)

bool puntoInternoPoligono(Vector3d& punto, const vector<Eigen::Vector3d>& poligono, Vector3d& normale) {
    unsigned int numVertices = poligono.size();
    double angoloTotale = 0.0;

    // Trasformazione di ogni punto del poligono e del punto da verificare
    vector<Eigen::Vector3d> poligonoTraslato;
    Matrix3d Q;
    Q << normale[1], -((normale[0] * normale[1]) / sqrt(std::pow(normale[0], 2) + pow(normale[1], 2))), normale[0],
        -normale[0], -((normale[1] * normale[2]) / sqrt(std::pow(normale[0], 2) + pow(normale[1], 2))), normale[1],
        0, sqrt(pow(normale[0], 2) + pow(normale[1], 2)), normale[2];

    for (const auto& vertice : poligono) {
        poligonoTraslato.push_back(Q * vertice);
    }
    Vector3d puntoTraslato = Q * punto;

    for (unsigned int i = 0; i < numVertices; ++i) {
        Vector3d v1 = poligonoTraslato[i] - puntoTraslato;
        Vector3d v2 = poligonoTraslato[(i + 1) % numVertices] - puntoTraslato;

        double angolo = std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
        angoloTotale += angolo;
    }

    // Se la somma degli angoli è 2*PI, il punto è all'interno del poligono.
    if (fabs(angoloTotale - 2 * M_PI) < 1e-9 || fabs(angoloTotale - M_PI) < 1e-9) {
        return true;
    }
    return false;
}


bool puntoInSegmento(Vector3d& p1, Vector3d& p2, Vector3d& p3)
{
    Vector3d prodotto_Vettoriale = (p2-p1).cross(p3-p1);
    // Controllo se punti collineari
    if (prodotto_Vettoriale[0] == 0 && prodotto_Vettoriale[1] == 0 && prodotto_Vettoriale[2] == 0)
    {
        return ((p3[0] >= min(p1[0], p2[0]) && p3[0] <= max(p1[0], p2[0]) &&
                 p3[1] >= min(p1[1], p2[1]) && p3[1] <= max(p1[1], p2[1]) &&
                 p3[2] >= min(p1[2], p2[2]) && p3[2] <= max(p1[2], p2[2])));
    }
    return false;
}


void calcolaTipologiaTracce(Struttura_DFN& DFN)
{
    for (unsigned int k = 0; k < DFN.numFratture; k++)
    {
        for (unsigned int i = 0; i < DFN.numTracce; i++)
        {
            if (k == DFN.Id_Fratture_Intersecanti[i][0] || k == DFN.Id_Fratture_Intersecanti[i][1])
            {
                Vector2i tipologia;
                for (unsigned int j = 0; j < 2; j++)
                {
                    Vector3d p_traccia = DFN.coordinateTraccia[i][j];
                    bool risposta = true;
                    for (unsigned int l = 0; l < DFN.numVertici[k]; l++)
                    {
                        Vector3d p1 = DFN.coordinateVertici[k][l];
                        Vector3d p2 = DFN.coordinateVertici[k][(l+1)%DFN.numVertici[k]];
                        if (puntoInSegmento(p1, p2, p_traccia))
                        {
                            risposta = false;
                        }
                    }
                    tipologia[j] = risposta;
                }
                if (tipologia[0] == false && tipologia[1] == false)
                {
                    DFN.tipoTraccia[k].push_back({i,false});
                }
                else
                {
                    DFN.tipoTraccia[k].push_back({i,true});
                }
            }
        }
    }
}


vector<Vector2i> riordinaTracce(vector<double>& lunghezza, vector<Vector2i>& tipo)
{
    vector<pair<double, unsigned int>> coppieIdLunghezza;
    for (unsigned int i = 0; i < lunghezza.size(); i++)
    {
        coppieIdLunghezza.push_back({lunghezza[i], i});
    }
    SortLibrary::MergeSort(coppieIdLunghezza);

    vector<Vector2i> tipoRiordinatoPassante;
    vector<Vector2i> tipoRiordinatoNonPassante;
    for(unsigned int j = 0; j < coppieIdLunghezza.size(); j++)
    {
        for (unsigned int i = 0; i < tipo.size(); i++)
        {
            if (tipo[i][1] == false && tipo[i][0] == coppieIdLunghezza[j].second)
            {
                tipoRiordinatoPassante.push_back(tipo[i]);
                break;
            }
            else if (tipo[i][1] == true && tipo[i][0] == coppieIdLunghezza[j].second)
            {
                tipoRiordinatoNonPassante.push_back(tipo[i]);
                break;
            }
        }
    }

    vector<Vector2i> tipoRiordinato;
    for(unsigned int i = 0; i < tipoRiordinatoPassante.size(); i++)
    {tipoRiordinato.push_back(tipoRiordinatoPassante[i]);}
    for(unsigned int i = 0; i < tipoRiordinatoNonPassante.size(); i++)
    {tipoRiordinato.push_back(tipoRiordinatoNonPassante[i]);}
    return tipoRiordinato;
}

void calcolaLunghezzaTracce(Struttura_DFN& DFN)
{
    for (unsigned int i = 0; i < DFN.numTracce; i++)
    {
        double lunghezza = (DFN.coordinateTraccia[i][1] - DFN.coordinateTraccia[i][0]).norm();
        DFN.lunghezzaTraccia.push_back(lunghezza);
    }
    for (unsigned int i = 0; i < DFN.tipoTraccia.size(); i++)
    {
        DFN.tipoTraccia[i] = riordinaTracce(DFN.lunghezzaTraccia, DFN.tipoTraccia[i]);
    }

}

}

