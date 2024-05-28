#pragma once

#include <iostream>
#include "GeometryDFN.hpp"
#include "MergeSort.hpp"

using namespace std;

namespace DFN_Library{

//funzione che legge i valori del file e riempie corretamente la struttura
bool ImportFratture(const string &filename,
                Struttura_DFN& DFN);

//funzione che calcola il centroide di una frattura passata per il suo Id
Vector3d calcolaCentroide(vector<Vector3d>& vertici);

//funzione che controlla se due fratture sono vicine a sufficienza per cui c'è un'intersezione
bool possibiliTracce(vector<Vector3d>& poligono1, vector<Vector3d>& poligono2, double tolleranza);

// funzione che mi calcola la normale al poligono
Vector3d normalePoligono(vector<Vector3d>& poligono);

// funzione che mi calcoli le tracce


// funzione originale che lega tutte le altre funzioni
void calcolaTracce(Struttura_DFN& DFN);

// Funzione che mi veda se il punto di intersezione trovato è all'interno del poligono
bool puntoInternoPoligono(const Vector3d& punto, vector<Vector3d>& poligono);

// Funzione che controlla se il punto p1 è interno al triangolo formato dagli altri tre punti
bool puntointriangolo(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3, const Vector3d& p4);

// Funzione che controlli se tre punti sono collineari, e se il punto p3 è interno al segmento p1-p2
bool puntoInSegmento(Vector3d& p1, Vector3d& p2, Vector3d& p3);

// Funzione che calcoli il tipo di traccia, se passante false, se non passante true
void calcolaTipologiaTracce(Struttura_DFN& DFN);

// Funzione che mi calcoli la lunghezza delle tracce
void calcolaLunghezzaTracce(Struttura_DFN& DFN);

// Funzione che riordini le tracce
vector<Vector2i> riordinaTracce(vector<double>& lunghezza, vector<Vector2i>& tipo);

// Prima funzione di output
bool OutputTracce(Struttura_DFN& DFN, const string& fileOutput);

// Seconda funzione output
bool OutputFratture(Struttura_DFN& DFN, const string& fileOutput);
}
