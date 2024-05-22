#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "DFN.hpp"

using namespace std;

namespace FractureLibrary{

bool ImportFracture(const string& filename, DFN& dfn);
Vector4d CalculatePlane(const Fracture& f);
double CalculateDistance(const Vector3d point1, const Vector3d point2);
bool CalculateTraces(const Fracture &f1, const Fracture &f2);
}
