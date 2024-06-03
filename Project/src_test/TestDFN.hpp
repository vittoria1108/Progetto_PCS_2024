#ifndef __TESTDFN_H
#define __TESTDFN_H

#include <gtest/gtest.h>
#include "DFN.hpp"
#include "Eigen/Eigen"
#include <iostream>

using namespace Eigen;
using namespace std;

namespace FractureLibrary {
//****************************************************************
TEST(DFNTEST, TestComputeDistance){

    Vector3d point1 = {1, -2, 1};
    Vector3d point2 = {7, 0, -8};

    double distance = sqrt(CalculateSquareDistance(point1, point2));

    EXPECT_EQ(distance, 11);
}

//****************************************************************
TEST(DFNTEST, TestCalculatePlaneFourVertices){

    Fracture f;
    f.NumberVertices = 4;

    MatrixXd verticesCoordinate(3, 4);

    verticesCoordinate << 2, 4, 5, 1,
                          -3, -2, 2, 7,
                          8, 7, 6, 7;

    f.VerticesCoordinates = verticesCoordinate;

    Vector4d plane = f.CalculatePlane();
    plane /= plane.norm();

    Vector4d rightPlane = {3, 1, 7, 59};
    rightPlane /= rightPlane.norm();

    EXPECT_EQ(plane, rightPlane);
}

TEST(DFNTEST, TestCalculatePlaneFiveVertices){

    Fracture f;
    f.NumberVertices = 5;

    MatrixXd verticesCoordinate(3, 5);

    verticesCoordinate << 1, 0.81, -0.31, -0.81, 0.31,
                          0, 0.51, 0.82, -0.51, -0.82,
                          0, 0.29, 0.48, -0.29, -0.48;
    f.VerticesCoordinates = verticesCoordinate;

    f.VerticesCoordinates = verticesCoordinate;

    Vector4d plane = f.CalculatePlane();
    plane /= plane.norm();

    Vector4d rightPlane = {0, -0.29, 0.51, 0};
    rightPlane /= rightPlane.norm();

    for(unsigned int i = 0; i < 4; i++){
        EXPECT_NEAR(plane[i], rightPlane[i], 0.1);
    }
}

//****************************************************************
TEST(DFNTEST, TestCalculateR){

    Fracture f;
    f.NumberVertices = 3;

    MatrixXd verticesCoordinate(3, 3);

    verticesCoordinate << -2, -5, -2,
                          -3, -2, 2,
                          1, 2, 0;
    f.VerticesCoordinates = verticesCoordinate;

    Vector3d barycentre = f.VerticesCoordinates.rowwise().mean();
    f.Barycentre = barycentre;

    double r = f.CalculateR();

    EXPECT_NEAR(r, 3.32, 0.01);
}

//****************************************************************
TEST(DFNTEST, TestIntersectionTrue){

    double tol = 10 * numeric_limits<double>::epsilon();

    Fracture f;
    f.NumberVertices = 3;

    MatrixXd verticesCoordinate(3, 3);

    verticesCoordinate << 2, 1, 4,
                          4, 1, 1,
                          3, 2, 1;
    f.VerticesCoordinates = verticesCoordinate;

    Vector3d p_r = {1, 1, 2};
    Vector3d t_r = {-2, -1.5, 0};

    Vector2d beta;

    bool intersection = f.IntersectsLine(p_r, t_r, beta, tol);

    EXPECT_TRUE(intersection);
}

TEST(DFNTEST, TestIntersectionVerticeFalse){

    double tol = 10 * numeric_limits<double>::epsilon();

    Fracture f;
    f.NumberVertices = 4;

    MatrixXd verticesCoordinate(3, 4);

    verticesCoordinate << -2.34, -6.83, -5.5, -1.02,
                          1.16, 2.49, 6.98, 5.65,
                          0, 0, 0, 0;
    f.VerticesCoordinates = verticesCoordinate;

    Vector3d p_r = {-6.83, 2.49, 0};
    Vector3d t_r = {2.83, -2.49, 0};

    Vector2d beta;

    bool intersection = f.IntersectsLine(p_r, t_r, beta, tol);

    EXPECT_FALSE(intersection);
}

TEST(DFNTEST, TestIntersectionFalse){

    double tol = 10 * numeric_limits<double>::epsilon();

    Fracture f;
    f.NumberVertices = 5;

    MatrixXd verticesCoordinate(3, 5);

    verticesCoordinate << 1, 0.81, -0.31, -0.81, 0.31,
                          0, 0.51, 0.82, -0.51, -0.82,
                          0, 0.29, 0.48, -0.29, -0.48;
    f.VerticesCoordinates = verticesCoordinate;

    Vector3d p_r = {-1.89, 0, 0};
    Vector3d t_r = {0.33, 1.6, 0};

    Vector2d beta;

    bool intersection = f.IntersectsLine(p_r, t_r, beta, tol);

    EXPECT_FALSE(intersection);
}

//****************************************************************
TEST(DFNTEST, TestIntersectionLineFalse){

    double tol = 10 * numeric_limits<double>::epsilon();

    Vector4d plane1 = {3, 1, 7, 59};
    Vector4d plane2 = {3, 1, 7, 24};

    Vector3d p_r;
    Vector3d t_r;

    bool intersection = FindIntersectionLine(plane1, plane2, p_r, t_r, tol);

    EXPECT_FALSE(intersection);
}

TEST(DFNTEST, TestIntersectionLineTrue){

    double tol = 10 * numeric_limits<double>::epsilon();

    Vector4d plane1 = {3, 1, 7, 59};
    Vector4d plane2 = {-1, 7, 5, 17};

    Vector3d p_r;
    Vector3d t_r;

    bool intersection = FindIntersectionLine(plane1, plane2, p_r, t_r, tol);

    p_r /= p_r.norm();
    t_r /= t_r.norm();

    Vector3d rightP = {2, -3, 8};
    rightP /= rightP.norm();

    Vector3d rightT = {-2, -1, 1};
    rightT /= rightT.norm();

    EXPECT_TRUE(intersection);

    for(unsigned int i = 0; i < 3; i++){
        EXPECT_NEAR(p_r[i], rightP[i], 0.5);
        EXPECT_NEAR(t_r[i], rightT[i], 0.5);
    }
}

//****************************************************************
TEST(DFNTEST, TestReadDFNTrue){

    double tol = 10 * numeric_limits<double>::epsilon();
    const string fileName = "C:/Users/user/OneDrive/Desktop/Progetto_PCS_2024/Project/src_test/Test_DFN.txt";
    DFN dfn;

    bool success = ReadDFN(fileName, dfn, tol);

    MatrixXd f1Vertices(3, 4);

    f1Vertices << 2, 4, 3.49, 1,
                  -3, -2, 7.35, 7,
                  8, 7, 6, 7;

    MatrixXd f2Vertices(3, 3);

    f2Vertices << 2, 6, -2,
                  2, 2, 2,
                  4, 7, 7;

    MatrixXd f3Vertices(3, 3);

    f3Vertices << 2, 3.22, 2,
                  2, 2, -3,
                  7, 4.915, 4;

    EXPECT_EQ(dfn.NumberFractures, 3);

    for(unsigned int i = 0; i < 3; i++)
    {
        for(unsigned int j = 0; j < 4; j++)
        {
            EXPECT_NEAR(f1Vertices(i, j), dfn.Fractures[0].VerticesCoordinates(i, j), 0.1);
        }
    }

    for(unsigned int i = 0; i < 3; i++)
    {
        for(unsigned int j = 0; j < 3; j++)
        {
            EXPECT_NEAR(f2Vertices(i, j), dfn.Fractures[1].VerticesCoordinates(i, j), 0.1);
            EXPECT_NEAR(f3Vertices(i, j), dfn.Fractures[2].VerticesCoordinates(i, j), 0.1);
        }
    }

    EXPECT_EQ(dfn.NumberTraces, 2);

    EXPECT_NEAR(dfn.Traces[0].Length, 1.07, 0.1);
    EXPECT_NEAR(dfn.Traces[1].Length, 2.43, 0.1);

    EXPECT_EQ(dfn.Fractures[0].nTraces.size(), 1);
    EXPECT_EQ(dfn.Fractures[0].pTraces.size(), 0);

    EXPECT_EQ(dfn.Fractures[1].nTraces.size(), 1);
    EXPECT_EQ(dfn.Fractures[1].pTraces.size(), 1);

    EXPECT_EQ(dfn.Fractures[2].nTraces.size(), 0);
    EXPECT_EQ(dfn.Fractures[2].pTraces.size(), 1);

    EXPECT_TRUE(success);
}

TEST(DFNTEST, TestReadDFNFalse){

    double tol = 10 * numeric_limits<double>::epsilon();
    const string fileName = "DFN_NonEsistente";
    DFN dfn;

    bool success = ReadDFN(fileName, dfn, tol);

    EXPECT_FALSE(success);
}

//****************************************************************
TEST(DFNTEST, TestImportFracture){

    const string fileName = "C:/Users/user/OneDrive/Desktop/Progetto_PCS_2024/Project/src_test/Test_Fracture.txt";
    DFN dfn;

    bool success = ImportFracture(fileName, dfn);

    MatrixXd vertices(3, 5);

    vertices << 2, 4, 6.5, 0.97, -0.6,
                2, 2, 3.8, 5.73, 3.9,
                0, 0, 1, 2.07, 1.05;

    Vector3d barycentre = vertices.rowwise().mean();

    EXPECT_TRUE(success);

    EXPECT_EQ(dfn.Fractures[0].Id, 0);
    EXPECT_EQ(dfn.Fractures[0].NumberVertices, 5);

    for(unsigned int i = 0; i < 3; i++)
    {
        for(unsigned int j = 0; j < 5; j++)
        {
            EXPECT_NEAR(dfn.Fractures[0].VerticesCoordinates(i, j), vertices(i, j), 0.1);
        }
    }

    for(unsigned int i = 0; i < 3; i++)
    {
        EXPECT_NEAR(dfn.Fractures[0].Barycentre[i], barycentre[i], 0.1);
    }
}

/* TODO: Test CalculateTraces -> 1) Passare due fratture che non si intersecano
 *                               2) Passare due fratture che generano traccia passante
 *                               3) Passare due fratture che generano traccia non passante
                                 4) (forse, in caso anche da modificare metodo) Passare fratture nello stesso piano con lato coincidente */

}

#endif
