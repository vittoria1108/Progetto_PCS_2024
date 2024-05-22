#ifndef DFN_HPP
#define DFN_HPP

#include <iostream>
#include <Eigen/Eigen>


using namespace std;
using namespace Eigen;


namespace FractureLibrary{

struct Trace{
    Vector2i Ids = {};
    MatrixXd EndPointsCoordinates = {};
    map<unsigned int, bool> Tips = {};
};

struct Fracture{ //struct con gli attributi di ogni frattura del dfn
    unsigned int Id = 0;
    unsigned int NumVertices = 0;
    MatrixXd VerticesCoordinates = {};
    Vector3d Barycentre = {};
    vector<Trace> Traces = {};
};

struct DFN{

    unsigned int NumberFractures = 0;
    vector<Fracture> Fractures = {}; // vettore che riferisce ad ognuna delle fratture ( con i relativi attributi)
};


}

#endif // DFN_HPP
