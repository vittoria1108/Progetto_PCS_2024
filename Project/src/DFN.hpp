#ifndef DFN_HPP
#define DFN_HPP

#include <iostream>
#include <Eigen/Eigen>


using namespace std;
using namespace Eigen;


namespace FractureLibrary{

struct Trace{

    unsigned int Id;
    Vector2i FracturesIds = {};
    MatrixXd EndpointsCoordinates = {};
    double Length;

};

struct Fracture{

    unsigned int Id = 0;
    unsigned int NumberVertices = 0;
    MatrixXd VerticesCoordinates = {};
    Vector3d Barycentre = {};
    map<unsigned int, bool> Tips = {};
    vector<Trace> nTraces = {};
    vector<Trace> pTraces = {};

};

struct DFN{

    unsigned int NumberFractures = 0;
    vector<Fracture> Fractures = {};

    unsigned int NumberTraces = 0;
    vector<Trace> Traces = {};

};

}

#endif // DFN_HPP
