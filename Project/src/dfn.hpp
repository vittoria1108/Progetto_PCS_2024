#ifndef DFN_HPP
#define DFN_HPP

#include <iostream>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;


namespace DFNLibrary{

struct Trace{

    unsigned int Id;
    Vector2i FracturesIds = {};
    MatrixXd EndpointsCoordinates;
    double Length;
};


struct Fracture{        //struct con gli attributi di ogni frattura del dfn

    unsigned int Id = 0;
    unsigned int NumberVertices = 0;
    MatrixXd VerticesCoordinates = {}; //mappa con chiave l'id della frattura e con dentro il vettore contenente le coordinate dei vertici
    Vector3d Barycentre = {};
    vector<Trace> npTraces = {};
    vector<Trace> pTraces = {};
    map<unsigned int, bool> Tips = {};
};


struct DFN{

    unsigned int NumberFractures = 0;
    vector<Fracture> Fractures;

    unsigned int NumberTraces = 0;
    vector<Trace> Traces = {};
};


}

#endif // DFN_HPP
