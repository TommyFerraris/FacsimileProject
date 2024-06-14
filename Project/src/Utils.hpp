#pragma once

#include "GeometryDFN.hpp"

using namespace std;

namespace DFN_Library{

//funzione che legge i valori del file e riempie la struttura del DFN con le informazioni sulle Fratture
bool ImportFratture(const string &filename,
                Struttura_DFN& DFN);

//funzione che calcola il centroide di una frattura passata per il suo Id
Vector3d calcolaCentroide(const vector<Vector3d>& vertici);

//funzione che controlla se due fratture sono vicine a sufficienza per cui c'è un'intersezione (ovvero una Traccia)
bool possibiliTracce(const vector<Vector3d>& poligono1, const vector<Vector3d>& poligono2, const double tolleranza);

// funzione che mi calcola la normale al poligono
Vector3d normalePoligono(const vector<Vector3d>& poligono);

// funzione che lega tutte le altre funzioni, ovvero funzione che calcola le Traccie del DFN
void calcolaTracce(Struttura_DFN& DFN);

// Funzione che controlla se il punto di intersezione trovato è all'interno del poligono
bool puntoInternoPoligono(const Vector3d& punto, vector<Vector3d>& poligono);

// Funzione che controlla se il punto p1 è interno al triangolo formato dagli altri tre punti
bool puntointriangolo(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3, const Vector3d& p4);

// Funzione che controlla se tre punti sono collineari, ovvero se il punto p3 è interno al segmento (p1-p2)
bool puntoInSegmento(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3);

// Funzione che calcola il tipo di traccia, se passante false, se non passante true
void calcolaTipologiaTracce(Struttura_DFN& DFN);

// Funzione che calcola la lunghezza delle tracce
void calcolaLunghezzaTracce(Struttura_DFN& DFN);

// Funzione che riordina le tracce
vector<Vector2i> riordinaTracce(const vector<double>& lunghezza, const vector<Vector2i>& tipo);

// Prima funzione di output (restituisce informazioni sulle Tracce e le Fratture che le generano)
bool OutputTracce(Struttura_DFN& DFN, const string& fileOutput);

// Seconda funzione output (restituisce informazioni sulla tipologia delle Tracce e di quante Tracce sono presenti in una Frattura)
bool OutputFratture(Struttura_DFN& DFN, const string& fileOutput);
}
