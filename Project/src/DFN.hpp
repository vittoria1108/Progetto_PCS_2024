#ifndef DFN_HPP
#define DFN_HPP

#include <iostream>
#include <Eigen/Eigen>


using namespace std;
using namespace Eigen;


namespace FractureLibrary{

struct Trace{

    Vector2i Ids = {};
    MatrixXd EndpointsCoordinates = {};
    map<unsigned int, bool> Tips = {};

};

struct Fracture{

    unsigned int Id = 0;
    unsigned int NumberVertices = 0;
    MatrixXd VerticesCoordinates = {};
    Vector3d Barycentre = {};
    vector<Trace> Traces = {};

};

struct DFN{

    unsigned int NumberFractures = 0;
    vector<Fracture> Fractures = {};

};

}

#endif // DFN_HPP
