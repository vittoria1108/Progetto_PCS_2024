#ifndef __POLYGONALMESH_H
#define __POLYGONALMESH_H

#include <iostream>
#include "Eigen/Eigen"
#include "DFN.hpp"

using namespace std;
using namespace Eigen;
using namespace FractureLibrary;

namespace PolygonalLibrary{


struct Cell0D {

    unsigned int Id = -1;
    Vector3d Coordinates = {};

};

struct Cell1D {

    unsigned int Id = -1;
    Vector2i Vertices = {};

    bool IsOld = false;
    unsigned int CuttedBy = -1;

    vector<unsigned int> NearCells2D = {};

};

struct Cell2D {

    unsigned int Id = -1;

    unsigned int NumberVertices = 0;
    vector<unsigned int> Vertices = {};

    unsigned int NumberEdges = 0;
    vector<unsigned int> Edges = {};

    bool IsOld = false;

};

struct PolygonalMesh {

    unsigned int NumberCell0D = 0;
    vector<Cell0D> Cells0D = {};

    unsigned int NumberCell1D = 0;
    vector<Cell1D> Cells1D = {};

    unsigned int NumberCell2D = 0;
    vector<Cell2D> Cells2D = {};

};

bool AlreadyExists(const Vector3d &coordinates,
                   const PolygonalMesh &PM,
                   unsigned int &id,
                   const double &tol);

bool IntersectionCellTrace(const Vector3d &p_s,
                           const Vector3d &t_s,
                           const Vector3d &p_r,
                           const Vector3d &t_r,
                           Vector3d &coordinates,
                           const double &tol);

bool IntersectionCellTrace(const Vector3d &p_s,
                           const Vector3d &t_s,
                           const Vector3d &p_r,
                           const Vector3d &t_r,
                           vector<double> &beta,
                           Vector3d &coordinates,
                           const double &tol);

void CreateFirstCell(PolygonalMesh &PM,
                     const Fracture &f,
                     unsigned int &idCell2D);

void ImportMesh(PolygonalMesh &PM,
                const Fracture &f,
                const double &tol);

}

#endif
