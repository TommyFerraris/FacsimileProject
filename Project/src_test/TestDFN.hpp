#pragma once

#include "GeometryDFN.hpp"
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include <gtest/gtest.h>
#include <vector>
#include "Eigen/Eigen"
#include "cmath"
#include "algorithm"

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

