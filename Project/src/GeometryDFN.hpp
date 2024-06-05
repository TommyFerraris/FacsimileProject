#pragma once

#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace DFN_Library {

struct Struttura_DFN{
    unsigned int numFratture = 0;
    vector<unsigned int> Id_Fratture = {}; /// dimensione 1 x numFratture
    vector<unsigned int> numVertici = {}; /// dimensione 1 x numFratture
    map <unsigned int, vector<Vector3d>> coordinateVertici = {}; /// marker Id_fratture, dimensione (3 x numVertici) x numFratture

    unsigned int numTracce = 0;
    vector<unsigned int> Id_Traccia = {}; /// dimensione 1 x numTracce
    vector<Vector2i> Id_Fratture_Intersecanti = {}; /// dimensione 2 x numTracce
    map <unsigned int, array<Vector3d, 2>> coordinateTraccia = {}; /// marker Id_traccia, dimensione (3 x 2) x numTracce
    vector<double> lunghezzaTraccia = {}; /// dimensione 1 x numTracce
    map <unsigned int, vector<Vector2i>> tipoTraccia = {}; /// dimensione 2 x (2 x numTracce)
    vector<unsigned int> Id_FrattureConTraccia = {}; /// dimensione 1 x dimensione TipoTraccia
};

}
