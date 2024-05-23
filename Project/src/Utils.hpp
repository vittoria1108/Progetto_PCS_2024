#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include "DFN.hpp"

using namespace std;

namespace FractureLibrary{

bool ImportFracture(const string &fileName,
                    DFN &dfn);

double CalculateDistance(const Vector3d point1,
                         const Vector3d point2);

Vector4d CalculatePlane(const Fracture &f);

bool CalculateTraces(Fracture &f1,
                     Fracture &f2,
                     unsigned int &id);

}




#endif // UTILS_HPP
