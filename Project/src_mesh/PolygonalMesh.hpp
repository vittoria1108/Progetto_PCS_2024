#ifndef __POLYGONALMESH_H
#define __POLYGONALMESH_H

#include <iostream>
#include <Eigen/Eigen>
#include <DFN.hpp>


using namespace std;
using namespace Eigen;
using namespace FractureLibrary;


namespace PolygonalLibrary{

struct Cell0D {

    unsigned int Id = -1;
    Eigen::Vector3d Coordinates = {};

};

struct Cell1D {

    unsigned int Id = -1;
    Eigen::Vector2i Vertices = {};

    bool IsOld = false;
    unsigned int CuttedBy = -1;

    std::list<unsigned int> NearCells2D = {};

};

struct Cell2D {

    unsigned int Id = -1;

    unsigned int NumberVertices = 0;
    std::vector<unsigned int> Vertices = {};

    unsigned int NumberEdges = 0;
    std::vector<unsigned int> Edges = {};

    bool IsOld = false;
    bool IsValid = true;

    double CalculateArea(const std::vector<Cell0D> &allCells0D)
    {
        double area = 0;

        unsigned int refPoint = Vertices[0];
        Eigen::Vector3d refPointCoordinates = allCells0D[refPoint].Coordinates;

        for(unsigned int i = 1; i < Vertices.size() - 1; i++)
        {
            unsigned int firstPoint = Vertices[i];
            unsigned int secondPoint = Vertices[i + 1];

            Eigen::Vector3d firstEdge = allCells0D[firstPoint].Coordinates - refPointCoordinates;
            Eigen::Vector3d secondEdge = allCells0D[secondPoint].Coordinates - refPointCoordinates;

            double areaTemp = 0.5 * (firstEdge.cross(secondEdge)).norm();
            area += areaTemp;
        }

        return area;
    }
};

struct PolygonalMesh {

    unsigned int NumberCell0D = 0;
    std::vector<Cell0D> Cells0D = {};

    unsigned int NumberCell1D = 0;
    std::vector<Cell1D> Cells1D = {};

    unsigned int NumberCell2D = 0;
    std::vector<Cell2D> Cells2D = {};

};

bool AlreadyExists(const Eigen::Vector3d &coordinates,
                   const PolygonalMesh &PM,
                   unsigned int &id,
                   const double &tol);

bool IntersectionCellTrace(const Eigen::Vector3d &p_s,
                           const Eigen::Vector3d &t_s,
                           const Eigen::Vector3d &p_r,
                           const Eigen::Vector3d &t_r,
                           Eigen::Vector3d &coordinates,
                           const double &tol);

bool CellContainsTrace(const PolygonalMesh &PM,
                       const Cell2D *cell,
                       const Eigen::Vector3d &p_r,
                       const Eigen::Vector3d &t_r,
                       Eigen::Vector2d &beta,
                       const double &tol);

void CreateFirstCell(PolygonalMesh &PM,
                     const FractureLibrary::Fracture &f,
                     unsigned int &idCell2D);

void CreateNewCellsByTrace(PolygonalMesh &PM,
                           unsigned int &idCell1D,
                           unsigned int &idCell2D,
                           const std::vector<Cell0D> &vertices,
                           const unsigned int &indexNewEdge,
                           const std::vector<unsigned int> &firstCell2DVertices,
                           const std::vector<unsigned int> &secondCell2DVertices,
                           std::vector<unsigned int> &firstCell2DEdges,
                           std::vector<unsigned int> &secondCell2DEdges,
                           const double &tol);

void CreateLastCell(PolygonalMesh &PM,
                    unsigned int &idCell2D,
                    const Cell2D *cell2D,
                    const vector<Cell1D> &newEdges);

void CreateNewCells(PolygonalMesh &PM,
                    const FractureLibrary::Trace &t,
                    unsigned int &idCell0D,
                    unsigned int &idCell1D,
                    unsigned int &idCell2D,
                    const bool &pass,
                    const double &tol);

void GenerateMesh(PolygonalMesh &PM,
                const FractureLibrary::Fracture &f,
                const double &tol);

}

#endif
