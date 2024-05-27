#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include "dfn.hpp"

using namespace std;

namespace DFNLibrary{

bool CompareTraces(const Trace &t1,
                   const Trace &t2);


// funzione per importare il file dei dati
bool ImportFracture(const string &fileName,
                DFN& dfn);


// funzione che calcola il piano
Vector4d CalculatePlane(const Fracture& fracture);


// calcola la distanza fra due punti
double CalculateDistance(const Vector3d point1,
                         const Vector3d point2);


// calcola i raggi delle fratture
double CalculateR(const Fracture &f);


// calcola l'intersezione fra i piani
bool FindIntersectionLine(const Vector4d &plane1,
                          const Vector4d &plane2,
                          Vector3d &p_r,
                          Vector3d &t_r,
                          const double &tol);


// calcola il segmento di intersezione delle fratture
bool IntersectionFractureLine(const Fracture &f,
                              const Vector3d &p_r,
                              const Vector3d &t_r,
                              Vector2d &beta,
                              const double &tol);


// calcola le tracce
void CalculateTraces(DFN &dfn,
                     Fracture &f1,
                     Fracture &f2,
                     unsigned int &id,
                     const double &tol);


void WriteOutputFiles(const string &outputTracesFile,
                      const string &outputTipsFile,
                      const DFN &dfn);


bool ReadDFN(const string &fileName,
             DFN &dfn,
             const double &tol);

}

#endif // UTILS_HPP
