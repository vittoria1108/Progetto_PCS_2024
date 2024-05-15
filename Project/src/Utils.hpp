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

bool CalculateTraces(const Fracture &f1,
                     const Fracture &f2);

}




#endif // UTILS_HPP
