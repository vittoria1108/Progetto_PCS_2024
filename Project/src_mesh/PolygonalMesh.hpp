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

    unsigned int Id = 0;
    Vector3d Coordinates = {};

    bool IsBetween(const Cell0D &cell1, const Cell0D &cell2, const double &tol)
    {
        Vector3d firstRow = cell1.Coordinates - Coordinates;
        Vector3d secondRow = cell2.Coordinates - Coordinates;

        MatrixXd matrix(3, 3);

        matrix << firstRow[0], firstRow[1], firstRow[2],
                  secondRow[0], secondRow[1], secondRow[2],
                  1, 1, 1;

        if(abs(matrix.determinant()) < tol)
            return true;

        return false;
    }
};

struct Cell1D {

    unsigned int Id = 0;
    Vector2i Vertices = {};
    bool IsOld = false;

    vector<unsigned int> NearCells2D = {};

    bool IntersectsLine(const Vector3d &p_s, const Vector3d &t_s, const Vector3d &p_r, const Vector3d &t_r, Vector3d &coordinates, const double &tol)
    {
        Vector3d prod = t_s.cross(t_r);

        if (prod.norm() > tol)
        {
            double alpha = (((p_r - p_s).cross(t_r)).dot(prod)) / (prod.dot(prod));

            if(alpha >= 0 && alpha < 1)
            {
                coordinates = p_s + alpha * t_s;
                return true;
            }
        }

        return false;
    }

};

struct Cell2D {

    unsigned int Id = {};

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

void ImportMesh(PolygonalMesh &PM, const Fracture &f, const double &tol);

}

#endif
