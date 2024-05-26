#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include "DFN.hpp"

using namespace std;

namespace FractureLibrary{

bool CompareTraces(const Trace &t1,
                   const Trace &t2);

bool ImportFracture(const string& filename,
                    DFN& dfn,
                    const double &tol);

Vector4d CalculatePlane(const Fracture& f);

double CalculateDistance(const Vector3d point1,
                         const Vector3d point2);

double CalculateR(const Fracture &f);

bool FindIntersectionLine(const Vector4d &plane1,
                          const Vector4d &plane2,
                          Vector3d &p_r,
                          Vector3d &t_r,
                          const double &tol);

bool IntersectionFractureLine(const Fracture &f,
                              const Vector3d &p_r,
                              const Vector3d &t_r,
                              Vector2d &beta,
                              const double &tol);

void CalculateTraces(DFN &dfn,
                     Fracture &f1,
                     Fracture &f2,
                     unsigned int &id,
                     const double &tol);

void WriteOutputFiles(const string &outputTracesFile,
                      const string &outputTipsFile,
                      const DFN &dfn);
}
