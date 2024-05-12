#ifndef DFN_HPP
#define DFN_HPP

#include <iostream>
#include <Eigen/Eigen>


using namespace std;
using namespace Eigen;


namespace FractureLibrary{

struct DFN{

    unsigned int NumberFractures = 0;
    vector<unsigned int> FracturesId = {};
    vector<Vector3d> VerticesCoordinates = {};
    unsigned int NumberVertices = 0;

};

}

#endif // DFN_HPP
