#pragma once

#include "GeometryDFN.hpp"
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include <gtest/gtest.h>
#include <vector>
#include "Eigen/Eigen"
#include "cmath"
#include "algorithm"
#include "MergeSort.hpp"

using namespace std;
using namespace DFN_Library;
using namespace DFN_PolygonalLibrary;

// Test sulla funzione che calcola il centroide di un poligono di n lati
TEST(CalcolaCentroideTest, poligono2D) {
    vector<Vector3d> poligono1 = {
        Vector3d(1, 0, 0),
        Vector3d(0, 1, 0),
        Vector3d(0, 0, 0),
        Vector3d(1, 1, 0)
    };
    Vector3d risultato1(0.5, 0.5, 0);
    ASSERT_EQ(calcolaCentroide(poligono1), risultato1);
}
TEST(CalcolaCentroideTest, poligono3D) {
    vector<Vector3d> poligono2 = {
        Vector3d(2, 3, -4),
        Vector3d(-2, -3, 4),
        Vector3d(2, -3, 4),
        Vector3d(-2, 3, -4)
    };
    Vector3d risultato2(0, 0, 0);
    ASSERT_EQ(calcolaCentroide(poligono2), risultato2);
}
TEST(CalcolaCentroideTest, poligono6Vertici) {
    vector<Vector3d> poligono3 = {
        Vector3d(1, 2, 3),
        Vector3d(4, 5, 6),
        Vector3d(7, 8, 9),
        Vector3d(10, 11, 12),
        Vector3d(13, 14, 15),
        Vector3d(16, 17, 18)
    };
    Vector3d risultato3(8.5, 9.5, 10.5);
    ASSERT_EQ(calcolaCentroide(poligono3), risultato3);
}


TEST(possibiliTracceTest, Poligonivicini){
    vector<Vector3d> poligono1 = {
        Vector3d(1, 0, 0),
        Vector3d(0, 1, 0),
        Vector3d(0, 0, 0),
        Vector3d(1, 1, 0)
    };
    vector<Vector3d> poligono2 = {
        Vector3d(0.7, -0.2, 0.3),
        Vector3d(0, 0, 0),
        Vector3d(-0.1, 0.5, 0.7),
        Vector3d(0, 0.2, 0.9)
    };
    const double tolleranza1 = 0.5;
    ASSERT_TRUE(possibiliTracce(poligono1, poligono2, tolleranza1));
}
TEST(possibiliTracceTest, Poligonilontani){
    vector<Vector3d> poligono1 = {
        Vector3d(1, 0, 0),
        Vector3d(0, 1, 0),
        Vector3d(0, 0, 0),
        Vector3d(1, 1, 0)
    };
    vector<Vector3d> poligono2 = {
        Vector3d(2.7, 2.2, 3.3),
        Vector3d(2, 2, 4),
        Vector3d(2.1, 2.5, 3.7),
        Vector3d(2, 2.2, 3.9),
        Vector3d(2.3, 2, 3.1),
    };
    const double tolleranza1 = 0.5;
    ASSERT_FALSE(possibiliTracce(poligono1, poligono2, tolleranza1));
}


TEST(NormalePoligonoTest, PoligonoPianoXY) {
    vector<Vector3d> poligono = {
        Vector3d(0, 0, 0),
        Vector3d(1, 0, 0),
        Vector3d(1, 1, 0)
    };
    Vector3d normale(0, 0, -1);
    ASSERT_EQ(normalePoligono(poligono), normale);
}
TEST(NormalePoligonoTest, PoligonoPianoXZ) {
    std::vector<Vector3d> poligono = {
        Vector3d(0, 0, 0),
        Vector3d(1, 0, 0),
        Vector3d(1, 0, 1)
    };
    Vector3d normale(0, 1, 0);
    ASSERT_EQ(normalePoligono(poligono), normale);
}

//Test per controllo della funzione che controlla se un punto è interno o meno ad un poligono
TEST(PuntoInternoPoligonoTest, PoligonoPianoXY){
    vector<Vector3d> poligono = {
        Vector3d(0, 0, 0),
        Vector3d(1, 0, 0),
        Vector3d(1, 1, 0),
        Vector3d(0, 1, 0)
    };
    Vector3d punto(0.5, 0.5, 0);
    ASSERT_TRUE(puntoInternoPoligono(punto, poligono));
}

TEST(PuntoInternoPoligonoTest, PoligonoPianoXZ){
    vector<Vector3d> poligono = {
        Vector3d(0, 0, 0),
        Vector3d(0, 0, 1),
        Vector3d(1, 0, 1),
        Vector3d(1, 0, 0)
    };
    Vector3d punto(0.5, 0, 0.5);
    ASSERT_TRUE(puntoInternoPoligono(punto, poligono));
}

TEST(PuntoInternoPoligonoTest, PoligonoPianoYZ){
    vector<Vector3d> poligono = {
        Vector3d(0, 0, 0),
        Vector3d(0, 1, 0),
        Vector3d(0, 1, 1),
        Vector3d(0, 0, 1)
    };
    Vector3d punto(0, 0.5, 0.5);
    ASSERT_TRUE(puntoInternoPoligono(punto, poligono));
}

//Test per controllare il punto p1 appartiane al triangolo formato da altri tre punti
TEST(PuntoInternoTriangoloTest, TriangoloXY){
    Vector3d p1(0.3, 0.3, 0);
    Vector3d p2(0, 0, 0);
    Vector3d p3(1, 0, 0);
    Vector3d p4(0, 1, 0);
    ASSERT_TRUE(puntointriangolo(p1, p2, p3, p4));
}
TEST(PuntoInternoTriangoloTest, TriangoloXZ){
    Vector3d p1(0.3, 0, 0.3);
    Vector3d p2(0, 0, 0);
    Vector3d p3(0, 0, 1);
    Vector3d p4(1, 0, 0);
    ASSERT_TRUE(puntointriangolo(p1, p2, p3, p4));
}
TEST(PuntoInternoTriangoloTest, TriangoloYZ){
    Vector3d p1(0, 0.3, 0.3);
    Vector3d p2(0, 0, 0);
    Vector3d p3(0, 1, 0);
    Vector3d p4(0, 0, 1);
    ASSERT_TRUE(puntointriangolo(p1, p2, p3, p4));
}

//Test sulla funzione che controlla se tre punti siano collineari, ovvero controlla che il punto p3 appartenga al segmento generato dai punti p1 e p2
TEST(PuntoInSegmentoTest, TestTrue){
    Vector3d p1(0, 0, 0);
    Vector3d p2(1, 1, 1);
    Vector3d p3(0.5, 0.5, 0.5);
    ASSERT_TRUE(puntoInSegmento(p1, p2, p3));
}
TEST(PuntoInSegmentoTest, TestFalse){
    Vector3d p1(0, 0, 0);
    Vector3d p2(1, 1, 1);
    Vector3d p3(0, 0, 1);
    ASSERT_FALSE(puntoInSegmento(p1, p2, p3));
}

// Test per controllare il riordinamento delle tracce
TEST(RiordinaTracceTest, TestBase){
    vector<double> lunghezza = {5.65, 10.0, 7.91};
    vector<Vector2i> tipo = {
        Vector2i{0, true},
        Vector2i{1, true},
        Vector2i{2, false}
    };

    // Ordine atteso: passanti ordinati per lunghezza + non passanti ordinati per lunghezza
    vector<Vector2i> vettoreOrdinato = {
        Vector2i{2, false}, // 10.0
        Vector2i{1, true},  // 7.91
        Vector2i{0, true}   // 5.65
    };
    vector<Vector2i> funzione = riordinaTracce(lunghezza, tipo);
    ASSERT_EQ(funzione.size(), vettoreOrdinato.size());
    for (size_t i = 0; i < funzione.size(); ++i) {
        EXPECT_EQ(funzione[i], vettoreOrdinato[i]) << "Differ at index " << i;
    }
}

TEST(RiordinaTracceTest, TestUgualeLunghezzaTipoDiverso){
    vector<double> lunghezza = {10.0, 10.0};
    vector<Vector2i> tipo = {
        Vector2i{0, true},
        Vector2i{1, false}
    };

    // Ordine atteso: passanti ordinati per lunghezza + non passanti ordinati per lunghezza
    vector<Vector2i> vettoreOrdinato = {
        Vector2i{1,false}, // 10.0 tipo passante
        Vector2i{0, true},  // 10.0 tipo non passante
    };
    vector<Vector2i> funzione = riordinaTracce(lunghezza, tipo);
    ASSERT_EQ(funzione.size(), vettoreOrdinato.size());
    for (size_t i = 0; i < funzione.size(); ++i) {
        EXPECT_EQ(funzione[i], vettoreOrdinato[i]) << "Differ at index " << i;
    }
}

// Test per calcolare lunghezza tracce
TEST(CalcolaLunghezzaTracceTest, TestBase) {
    Struttura_DFN DFN;
    DFN.numTracce = 2;
    DFN.coordinateTraccia[0] = {Vector3d{0, 0, 0}, Vector3d{1, 0, 0}};
    DFN.coordinateTraccia[1] = {Vector3d{1, 1, 1}, Vector3d{2, 0, 0}};
    DFN.tipoTraccia[2] = {Vector2i{0, false}, Vector2i{1, true}};

    calcolaLunghezzaTracce(DFN);

    ASSERT_EQ(DFN.lunghezzaTraccia.size(), 2);
    EXPECT_DOUBLE_EQ(DFN.lunghezzaTraccia[0], 1.0);
    EXPECT_DOUBLE_EQ(DFN.lunghezzaTraccia[1], sqrt(3.0));
}

TEST(CalcolaLunghezzaTracceTest, TestLunghezzaNulla) {
    Struttura_DFN DFN;
    DFN.numTracce = 2;
    DFN.coordinateTraccia[0] = {Vector3d{0, 0, 0}, Vector3d{0, 0, 0}};
    DFN.coordinateTraccia[1] = {Vector3d{1, 1, 1}, Vector3d{4, 5, 1}};
    DFN.tipoTraccia[2] = {Vector2i{0, false}, Vector2i{1, true}};

    calcolaLunghezzaTracce(DFN);

    ASSERT_EQ(DFN.lunghezzaTraccia.size(), 2);
    EXPECT_DOUBLE_EQ(DFN.lunghezzaTraccia[0], 0.0);
    EXPECT_DOUBLE_EQ(DFN.lunghezzaTraccia[1], 5.0);
}

// Funzione che testi se calcolo in maniera correta la tipologia delle tracce e la salvo
TEST(CalcolaTipologiaTracceTest, FunzioneBase) {
    Struttura_DFN DFN;
    DFN.numFratture = 2;
    DFN.numTracce = 3;
    DFN.numVertici = {3, 4};

    // Coordinate delle tracce
    DFN.coordinateTraccia[0] ={Vector3d(1, 1, 0), Vector3d(2, 2, 0)};
    DFN.coordinateTraccia[1] ={Vector3d(3, 3, 0), Vector3d(4, 4, 0)};
    DFN.coordinateTraccia[2] ={Vector3d(5, 5, 0), Vector3d(6, 6, 0)};

    // Coordinate dei vertici delle fratture
    DFN.coordinateVertici[0] = {Vector3d(0, 0, 0), Vector3d(2, 0, 0), Vector3d(1, 1, 0)};
    DFN.coordinateVertici[1] = {Vector3d(3, 3, 0), Vector3d(5, 3, 0), Vector3d(5, 5, 0), Vector3d(3, 5, 0)};

    // Id delle fratture intersecanti
    DFN.Id_Fratture_Intersecanti = {
        {0, 1},
        {0, 1},
        {1, 0}
    };

    // Chiamata della funzione da testare
    calcolaTipologiaTracce(DFN);

    // Fratture con traccia attese
    vector<int> Id_FrattureConTraccia_Expected = {0, 1};

    // Verifica
    ASSERT_EQ(DFN.Id_FrattureConTraccia.size(), Id_FrattureConTraccia_Expected.size());
    for (size_t i = 0; i < DFN.Id_FrattureConTraccia.size(); ++i) {
        EXPECT_EQ(DFN.Id_FrattureConTraccia[i], Id_FrattureConTraccia_Expected[i]);
    }

    // Verifica delle tipologie delle tracce
    map<unsigned int, vector<Vector2i>> tipoTraccia_Expected;
    tipoTraccia_Expected[0] = {{0, true}, {1, true}, {2, true}};
    tipoTraccia_Expected[1] = {{0, true}, {1, true}, {2, true}};

    ASSERT_EQ(DFN.tipoTraccia.size(), tipoTraccia_Expected.size());
    for (size_t i = 0; i < DFN.tipoTraccia.size(); ++i) {
        ASSERT_EQ(DFN.tipoTraccia[i].size(), tipoTraccia_Expected[i].size());
        for (size_t j = 0; j < DFN.tipoTraccia[i].size(); ++j) {
            EXPECT_EQ(DFN.tipoTraccia[i][j], tipoTraccia_Expected[i][j]) << "Differ at index " << i << ", " << j;
        }
    }
}

TEST(CalcolaTipologiaTracceTest, FunzioneConPassante) {
    Struttura_DFN DFN;
    DFN.numFratture = 2;
    DFN.numTracce = 1;
    DFN.numVertici = {4, 4};

    // Coordinate delle tracce
    DFN.coordinateTraccia[0] ={Vector3d(0.8, 0, 0), Vector3d(0.8, 1, 0)};

    // Coordinate dei vertici delle fratture
    DFN.coordinateVertici[0] = {Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(1, 1, 0), Vector3d(0, 1, 0)};
    DFN.coordinateVertici[1] = {Vector3d(0.8, 0, -1), Vector3d(0.8, 0, 0.3), Vector3d(0.8, 1, 0.3), Vector3d(0.8, 1, -1)};

    // Id delle fratture intersecanti
    DFN.Id_Fratture_Intersecanti = {{0, 1}};

    // Chiamata della funzione da testare
    calcolaTipologiaTracce(DFN);

    // Fratture con traccia attese
    vector<int> Id_FrattureConTraccia_Expected = {0, 1};

    // Verifica
    ASSERT_EQ(DFN.Id_FrattureConTraccia.size(), Id_FrattureConTraccia_Expected.size());
    for (size_t i = 0; i < DFN.Id_FrattureConTraccia.size(); ++i) {
        EXPECT_EQ(DFN.Id_FrattureConTraccia[i], Id_FrattureConTraccia_Expected[i]);
    }

    // Verifica delle tipologie delle tracce
    map<unsigned int, vector<Vector2i>> tipoTraccia_Expected;
    tipoTraccia_Expected[0] = {{0, false}};
    tipoTraccia_Expected[1] = {{0, false}};

    ASSERT_EQ(DFN.tipoTraccia.size(), tipoTraccia_Expected.size());
    for (size_t i = 0; i < DFN.tipoTraccia.size(); ++i) {
        ASSERT_EQ(DFN.tipoTraccia[i].size(), tipoTraccia_Expected[i].size());
        for (size_t j = 0; j < DFN.tipoTraccia[i].size(); ++j) {
            EXPECT_EQ(DFN.tipoTraccia[i][j], tipoTraccia_Expected[i][j]) << "Differ at index " << i << ", " << j;
        }
    }
}

