#ifndef __TESTDFN_H
#define __TESTDFN_H

#include <gtest/gtest.h>
#include "DFN.hpp"
#include "Eigen/Eigen"
#include <iostream>

using namespace Eigen;
using namespace std;

namespace FractureLibrary {
//********************************
TEST(FRACTURETEST, TestComputeDistance){

    Vector3d point1 = {1, 0, 0};
    Vector3d point2 = {7, 0, 0};

    double distance = CalculateSquareDistance(point1, point2);

    EXPECT_EQ(distance, 36);
}


}

#endif
