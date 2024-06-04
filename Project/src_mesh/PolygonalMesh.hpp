#pragma once
#include <iostream>
#include <Eigen/Eigen>
#include "DFN.hpp"

using namespace std;
using namespace Eigen;
using namespace DFNLibrary;

namespace PolygonalLibrary{

struct PolygonalMesh {

    unsigned int NumberCell0D = 0;
    vector<unsigned int> IdCell0D = {};
    vector<Vector3d> CoordinatesCell0D = {};

    unsigned int NumberCell1D = 0;
    vector<unsigned int> IdCell1D = {};
    vector<bool> Cell1DOld = {};
    vector<Vector2i> VerticesCell1D = {};

    unsigned int NumberCell2D = 0;
    vector<unsigned int> IdCell2D = {};
    vector<bool> Cell2DOld = {};
    vector<vector<unsigned int>> VerticesCell2D = {};
    vector<vector<unsigned int>> EdgesCell2D = {};
};

void ImportMesh(PolygonalMesh &PM,const Fracture &f);

}


