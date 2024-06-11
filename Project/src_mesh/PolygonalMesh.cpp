#include <iostream>
#include <Eigen/Eigen>
#include "DFN.hpp"
#include "PolygonalMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace FractureLibrary;

namespace PolygonalLibrary{

bool AlreadyExists(const Vector3d &coordinates,
                   const PolygonalMesh &PM,
                   unsigned int &id,
                   const double &tol)
{
    for(Cell0D cell : PM.Cells0D)
    {
        if(abs(cell.Coordinates[0] - coordinates[0]) < tol &&
            abs(cell.Coordinates[1] - coordinates[1]) < tol &&
            abs(cell.Coordinates[2] - coordinates[2]) < tol)
        {
            id = cell.Id;
            return true;
        }
    }

    return false;
}

bool IntersectionCellTrace(const Vector3d &p_s,
                           const Vector3d &t_s,
                           const Vector3d &p_r,
                           const Vector3d &t_r,
                           Vector3d &coordinates,
                           const double &tol)
{
    Vector3d prod = t_s.cross(t_r);

    if (prod.norm() > tol)
    {
        double alpha = (((p_r - p_s).cross(t_r)).dot(prod)) / (prod.dot(prod));

        if(alpha > 0 && alpha <= 1)
        {
            coordinates = p_s + (alpha * t_s);
            return true;
        }
    }

    return false;
}

bool IntersectionCellTrace(const Vector3d &p_s,
                           const Vector3d &t_s,
                           const Vector3d &p_r,
                           const Vector3d &t_r,
                           vector<double> &beta,
                           Vector3d &coordinates,
                           const double &tol)
{
    Vector3d prod = t_s.cross(t_r);

    if (prod.norm() > tol)
    {
        double alpha = (((p_r - p_s).cross(t_r)).dot(prod)) / (prod.dot(prod));

        if(alpha > 0 && alpha <= 1)
        {
            coordinates = p_s + (alpha * t_s);

            double betaTemp = (((p_s - p_r).cross(t_s)).dot(-prod)) / (prod.dot(prod));
            beta.push_back(betaTemp);

            return true;
        }
    }

    return false;
}

void CreateFirstCell(PolygonalMesh &PM,
                     const Fracture &f,
                     unsigned int &idCell2D)
{
    Cell2D fracture;
    fracture.Id = idCell2D++;
    fracture.NumberVertices = f.NumberVertices;
    fracture.NumberEdges = f.NumberVertices;

    for(unsigned int i = 0; i < f.NumberVertices; i++)
    {
        unsigned int next = i == f.NumberVertices - 1 ? 0 : i + 1;

        Cell0D cell0D;
        cell0D.Id = i;
        cell0D.Coordinates = f.VerticesCoordinates.col(i);
        fracture.Vertices.push_back(cell0D.Id);

        PM.NumberCell0D++;
        PM.Cells0D.push_back(cell0D);

        Cell1D cell1D;
        cell1D.Id = i;
        cell1D.Vertices = {i, next};
        cell1D.NearCells2D.push_back(fracture.Id);
        fracture.Edges.push_back(cell1D.Id);

        PM.NumberCell1D++;
        PM.Cells1D.push_back(cell1D);
    }

    PM.NumberCell2D++;
    PM.Cells2D.push_back(fracture);
}

void CreateNewCells(PolygonalMesh &PM,
                    const Trace &t,
                    unsigned int &idCell0D,
                    unsigned int &idCell1D,
                    unsigned int &idCell2D,
                    const bool &pass,
                    const double &tol)
{
    Vector3d p_r = t.EndpointsCoordinates.row(0);
    Vector3d t_r = t.EndpointsCoordinates.row(1) - p_r.transpose();

    unsigned int indexCell1D = 0;
    unsigned int indexCell2D = 0;

    while(PM.Cells1D[indexCell1D].IsOld)
    {
        indexCell1D++;
    }

    vector<unsigned int> cells2D = PM.Cells1D[indexCell1D].NearCells2D;

    while(cells2D.size() != 0)
    {
        indexCell2D = cells2D[0];

        if(PM.Cells2D[indexCell2D].IsOld)
        {
            cells2D.erase(cells2D.begin());
            continue;
        }

        vector<unsigned int> indexCell1DIntersected = {};

        Vector2i vertices;
        unsigned int indexVertices = 0;
        unsigned int indexNewVertice = 0;

        Cell2D cell2D = PM.Cells2D[indexCell2D];

        bool first = true;

        vector<unsigned int> firstCell2DVertices = {};
        vector<unsigned int> secondCell2DVertices = {};

        vector<unsigned int> firstCell2DEdges = {};
        vector<unsigned int> secondCell2DEdges = {};

        vector<Cell1D> newEdges = {};

        vector<double> beta;
        beta.reserve(2);

        for(unsigned int i : cell2D.Edges)
        {
            Cell1D cell1D = PM.Cells1D[i];

            Cell0D firstCell0D;
            Cell0D secondCell0D;

            for(unsigned int j = 0; j < cell2D.Vertices.size(); j++)
            {
                int next = j == cell2D.Vertices.size() - 1 ? 0 : j + 1;

                int firstVertice = cell2D.Vertices[j];

                auto it = find(cell1D.Vertices.begin(), cell1D.Vertices.end(), firstVertice);

                if (it == cell1D.Vertices.end())
                    continue;

                int secondVertice = cell2D.Vertices[next];

                it = find(cell1D.Vertices.begin(), cell1D.Vertices.end(), secondVertice);

                if (it == cell1D.Vertices.end())
                    continue;

                firstCell0D = PM.Cells0D[firstVertice];
                secondCell0D = PM.Cells0D[secondVertice];
                break;
            }

            Vector3d p_s = firstCell0D.Coordinates;
            Vector3d t_s = secondCell0D.Coordinates - p_s;

            Vector3d coordinates = {};

            if(first)
                firstCell2DVertices.push_back(firstCell0D.Id);
            else
                secondCell2DVertices.push_back(firstCell0D.Id);

            bool intersection;

            if(pass)
            {
                intersection = IntersectionCellTrace(p_s, t_s, p_r, t_r, coordinates, tol);
            }
            else
            {
                intersection = IntersectionCellTrace(p_s, t_s, p_r, t_r, beta, coordinates, tol);
            }

            if(intersection)
            {
                unsigned int id;
                Cell0D newCell0D;
                Cell1D firstCell1D;
                Cell1D secondCell1D;

                if(AlreadyExists(coordinates, PM, id, tol))
                {
                    newCell0D = PM.Cells0D[id];
                    vertices[indexVertices++] = newCell0D.Id;

                    Vector2i firstCellVertices = {firstCell0D.Id, newCell0D.Id};

                    bool found = false;

                    for(Cell1D &cell : PM.Cells1D)
                    {
                        if(found)
                        {
                            secondCell1D = cell; // lato successivo
                            break;
                        }

                        if(cell.Vertices == firstCellVertices)
                        {
                            firstCell1D = cell;
                            found = true;
                        }
                    }
                }
                else
                {
                    newCell0D.Id = idCell0D++;
                    newCell0D.Coordinates = coordinates;

                    PM.NumberCell0D++;
                    PM.Cells0D.push_back(newCell0D);

                    vertices[indexVertices++] = newCell0D.Id;

                    PM.Cells1D[i].IsOld = true;

                    firstCell1D.Id = idCell1D++;
                    firstCell1D.Vertices = {firstCell0D.Id, newCell0D.Id};

                    PM.NumberCell1D++;
                    PM.Cells1D.push_back(firstCell1D);

                    newEdges.push_back(firstCell1D);

                    secondCell1D.Id = idCell1D++;
                    secondCell1D.Vertices = {newCell0D.Id, secondCell0D.Id};

                    PM.NumberCell1D++;
                    PM.Cells1D.push_back(secondCell1D);

                    newEdges.push_back(secondCell1D);
                }

                PM.Cells1D[i].CuttedBy = newCell0D.Id;

                firstCell2DVertices.push_back(newCell0D.Id);
                secondCell2DVertices.push_back(newCell0D.Id);

                if(first)
                {
                    firstCell2DEdges.push_back(firstCell1D.Id);

                    indexNewVertice = firstCell2DEdges.size();
                    firstCell2DEdges.push_back(firstCell1D.Id);

                    secondCell2DEdges.push_back(secondCell1D.Id);
                }
                else
                {
                    firstCell2DEdges.push_back(secondCell1D.Id);
                    secondCell2DEdges.push_back(firstCell1D.Id);

                    indexCell1DIntersected.push_back(cell1D.Id);
                }

                first = first ? false : true;
            }
            else
            {
                if(first)
                    firstCell2DEdges.push_back(cell1D.Id);
                else
                    secondCell2DEdges.push_back(cell1D.Id);
            }
        }

        if(firstCell2DEdges != PM.Cells2D[indexCell2D].Edges)
        {
            Cell1D trace;
            trace.Id = idCell1D++;
            trace.Vertices = vertices;

            firstCell2DEdges[indexNewVertice] = trace.Id;
            secondCell2DEdges.push_back(trace.Id);

            Cell2D firstCell2D;
            firstCell2D.Id = idCell2D++;
            firstCell2D.NumberVertices = firstCell2DVertices.size();
            firstCell2D.Vertices = firstCell2DVertices;
            firstCell2D.NumberEdges = firstCell2DEdges.size();
            firstCell2D.Edges = firstCell2DEdges;

            PM.NumberCell2D++;
            PM.Cells2D.push_back(firstCell2D);

            Cell2D secondCell2D;
            secondCell2D.Id = idCell2D++;
            secondCell2D.NumberVertices = secondCell2DVertices.size();
            secondCell2D.Vertices = secondCell2DVertices;
            secondCell2D.NumberEdges = secondCell2DEdges.size();
            secondCell2D.Edges = secondCell2DEdges;

            PM.NumberCell2D++;
            PM.Cells2D.push_back(secondCell2D);

            for(unsigned int i : firstCell2DEdges)
            {
                if(i != trace.Id)
                    PM.Cells1D[i].NearCells2D.push_back(firstCell2D.Id);
            }

            for(unsigned int i : secondCell2DEdges)
            {
                if(i != trace.Id)
                    PM.Cells1D[i].NearCells2D.push_back(secondCell2D.Id);
            }

            trace.NearCells2D.push_back(firstCell2D.Id);
            trace.NearCells2D.push_back(secondCell2D.Id);

            PM.NumberCell1D++;
            PM.Cells1D.push_back(trace);

            PM.Cells2D[indexCell2D].IsOld = true;
            cells2D.erase(cells2D.begin());

            if(cells2D.size() == 0 && indexCell1DIntersected.size() != 0)
            {
                unsigned int index = indexCell1DIntersected[0] != indexCell1D ? indexCell1DIntersected[0] : indexCell1DIntersected[1];
                cells2D = PM.Cells1D[index].NearCells2D;
            }
        }
        else
        {
            if(cells2D.size() == 1)
                cells2D = PM.Cells1D[++indexCell1D].NearCells2D;
            else
                cells2D.erase(cells2D.begin());
        }

        if(!pass && (beta[0] >= 0 && beta[0] < 1) && beta[1] >= 1)
        {
            indexCell2D = cells2D[0];

            while(PM.Cells2D[indexCell2D].IsOld)
            {
                cells2D.erase(cells2D.begin());
                indexCell2D = cells2D[0];
            }

            cell2D = PM.Cells2D[indexCell2D];

            vector<unsigned int> newCell2DVertices = {};
            vector<unsigned int> newCell2DEdges = {};

            for(unsigned int i : cell2D.Edges)
            {
                Cell1D cell = PM.Cells1D[i];

                if(cell.IsOld)
                {
                    newCell2DVertices.push_back(cell.Vertices[0]);
                    newCell2DVertices.push_back(cell.CuttedBy);
                    newCell2DVertices.push_back(cell.Vertices[1]);

                    unsigned int firstEdge = 0;
                    unsigned int secondEdge = 0;

                    for(Cell1D newEdge : newEdges)
                    {
                        if(newEdge.Vertices[0] == cell.Vertices[0] && newEdge.Vertices[1] == cell.CuttedBy)
                        {
                            firstEdge = newEdge.Id;
                        }

                        if(newEdge.Vertices[0] == cell.CuttedBy && newEdge.Vertices[1] == cell.Vertices[1])
                        {
                            secondEdge = newEdge.Id;
                        }
                    }

                    newCell2DEdges.push_back(firstEdge);
                    newCell2DEdges.push_back(secondEdge);
                }
                else
                {
                    newCell2DVertices.push_back(cell.Vertices[0]);
                    newCell2DEdges.push_back(i);
                }

                PM.Cells2D[indexCell2D].IsOld = true;

                Cell2D newCell2D;
                newCell2D.Id = idCell2D++;
                newCell2D.NumberVertices = newCell2DVertices.size();
                newCell2D.Vertices = newCell2DVertices;
                newCell2D.NumberEdges = newCell2DEdges.size();
                newCell2D.Edges = newCell2DEdges;

                PM.NumberCell2D++;
                PM.Cells2D.push_back(newCell2D);
            }

            break;
        }
    }
}

void ImportMesh(PolygonalMesh &PM,
                const Fracture &f,
                const double &tol)
{

    unsigned int idCell2D = 0;

    unsigned int idCell0D = f.NumberVertices;
    unsigned int idCell1D = f.NumberVertices;

    // Memorizzo i dati della frattura come celle

    CreateFirstCell(PM, f, idCell2D);


    // Inizio a creare le celle nuove

    bool pass = true;

    for(const Trace &t : f.pTraces)
    {
        if(t.IsOnEdge.at(f.Id))
            continue;

        CreateNewCells(PM, t, idCell0D, idCell1D, idCell2D, pass, tol);
    }

    pass = false;

    for(const Trace &t : f.nTraces)
    {
        if(t.IsOnEdge.at(f.Id))
            continue;

        CreateNewCells(PM, t, idCell0D, idCell1D, idCell2D, pass, tol);
    }
}

}
