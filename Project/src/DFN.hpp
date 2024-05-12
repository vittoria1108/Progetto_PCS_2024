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
    map<int, vector<Vector3d>> VerticesCoordinates = {}; //mappa con chiave l'id della frattura e con dentro il vettore contenente le coordiante dei vertici
    vector<unsigned int> NumberVertices = {};

};

}

#endif // DFN_HPP
