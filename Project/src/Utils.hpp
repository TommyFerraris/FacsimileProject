#pragma once

#include <iostream>
#include "GeometryDFN.hpp"

using namespace std;

namespace DFN_Library{

//funzione che legge i valori del file e riempie corretamente la struttura
bool ImportFratture(const string &filename,
                Struttura_DFN& DFN);

//funzione che calcola il centroide di una frattura passata per il suo Id
Vector3d calcolaCentroide(Struttura_DFN& DFN);
}
