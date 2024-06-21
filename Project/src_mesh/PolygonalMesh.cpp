#include <iostream>
#include <Eigen/Eigen>
#include "DFN.hpp"
#include "PolygonalMesh.hpp"


using namespace std;
using namespace Eigen;
using namespace DFNLibrary;


namespace PolygonalLibrary{

bool AlreadyExists(const Eigen::Vector3d &coordinates,
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

bool IntersectionCellTrace(const Eigen::Vector3d &p_s,
                           const Eigen::Vector3d &t_s,
                           const Eigen::Vector3d &p_r,
                           const Eigen::Vector3d &t_r,
                           Eigen::Vector3d &coordinates,
                           const double &tol)
{
    Eigen::Vector3d prod = t_s.cross(t_r);

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

bool CellContainsTrace(const PolygonalMesh &PM,
                       const Cell2D *cell,
                       const Eigen::Vector3d &p_r,
                       const Eigen::Vector3d &t_r,
                       Eigen::Vector2d &beta,
                       const double &tol)
{
    unsigned int index = 0;

    bool intersection = false;

    for(unsigned int i : cell->Edges)
    {
        unsigned int firstCell0D = PM.Cells1D[i].Vertices[0];
        unsigned int secondCell0D = PM.Cells1D[i].Vertices[1];

        Eigen::Vector3d p_s = PM.Cells0D[firstCell0D].Coordinates;
        Eigen::Vector3d t_s = PM.Cells0D[secondCell0D].Coordinates - p_s;

        Eigen::Vector3d prod = t_s.cross(t_r);

        if (prod.norm() > tol)
        {
            double alpha = (((p_r - p_s).cross(t_r)).dot(prod)) / (prod.dot(prod));

            if(alpha > 0 && alpha <= 1)
            {
                intersection = true;

                double betaTemp = (((p_s - p_r).cross(t_s)).dot(-prod)) / (prod.dot(prod));
                beta[index++] = betaTemp;

                if(index == 2)
                    break;
            }
        }
    }

    if(!intersection)
        return false;

    if((beta[0] <= 0 && beta[1] <= 0) || (beta[0] >= 1 && beta[1] >= 1))
        return false;

    sort(beta.begin(), beta.end());
    return true;
}


void CreateFirstCell(PolygonalMesh &PM,
                     const DFNLibrary::Fracture &f,
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

void CreateNewCellsByTrace(PolygonalMesh &PM,
                           unsigned int &idCell1D,
                           unsigned int &idCell2D,
                           const std::vector<Cell0D> &vertices,
                           const unsigned int &indexNewEdge,
                           const std::vector<unsigned int> &firstCell2DVertices,
                           const std::vector<unsigned int> &secondCell2DVertices,
                           std::vector<unsigned int> &firstCell2DEdges,
                           std::vector<unsigned int> &secondCell2DEdges,
                           const double &tol)
{
    Cell1D trace;
    trace.Id = idCell1D++;

    Eigen::Vector2i traceVertices = {vertices[0].Id, vertices[1].Id};
    trace.Vertices = traceVertices;

    firstCell2DEdges[indexNewEdge] = trace.Id;
    secondCell2DEdges.push_back(trace.Id);

    Cell2D firstCell2D;
    firstCell2D.Id = idCell2D++;
    firstCell2D.NumberVertices = firstCell2DVertices.size();
    firstCell2D.Vertices = firstCell2DVertices;
    firstCell2D.NumberEdges = firstCell2DEdges.size();
    firstCell2D.Edges = firstCell2DEdges;

    double firstArea = firstCell2D.CalculateArea(PM.Cells0D);

    if(firstArea < tol*tol)
    {
        firstCell2D.IsValid = false;
    }

    Cell2D secondCell2D;
    secondCell2D.Id = idCell2D++;
    secondCell2D.NumberVertices = secondCell2DVertices.size();
    secondCell2D.Vertices = secondCell2DVertices;
    secondCell2D.NumberEdges = secondCell2DEdges.size();
    secondCell2D.Edges = secondCell2DEdges;

    double secondArea = secondCell2D.CalculateArea(PM.Cells0D);

    if(secondArea < tol*tol)
    {
        secondCell2D.IsValid = false;
    }

    PM.NumberCell2D++;
    PM.Cells2D.push_back(firstCell2D);

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
}

void CreateLastCell(PolygonalMesh &PM,
                    unsigned int &idCell2D,
                    const Cell2D *cell2D)
{

    Cell2D newCell2D;
    newCell2D.Id = idCell2D++;

    std::vector<unsigned int> newCell2DVertices = {};
    std::vector<unsigned int> newCell2DEdges = {};

    for(unsigned int i : cell2D->Edges)
    {
        Cell1D *cell = &PM.Cells1D[i];

        int firstVertice;
        int secondVertice;

        for(unsigned int j = 0; j < cell2D->Vertices.size(); j++)
        {
            int next = j == cell2D->Vertices.size() - 1 ? 0 : j + 1;

            firstVertice = cell2D->Vertices[j];

            auto it = find(cell->Vertices.begin(), cell->Vertices.end(), firstVertice);

            if (it == cell->Vertices.end())
                continue;

            secondVertice = cell2D->Vertices[next];

            it = find(cell->Vertices.begin(), cell->Vertices.end(), secondVertice);

            if(it == cell->Vertices.end())
                continue;

            break;
        }

        newCell2DVertices.push_back(firstVertice);

        if(cell->IsOld)
        {
            unsigned int firstEdgeId = cell->ReplacedBy[0];
            unsigned int secondEdgeId = cell->ReplacedBy[1];

            Cell1D *firstEdge = &PM.Cells1D[firstEdgeId];
            Cell1D *secondEdge = &PM.Cells1D[secondEdgeId];

            auto it = find(firstEdge->Vertices.begin(), firstEdge->Vertices.end(), firstVertice);

            if(it == firstEdge->Vertices.end())
            {
                firstEdgeId = cell->ReplacedBy[1];
                secondEdgeId = cell->ReplacedBy[0];

                firstEdge = &PM.Cells1D[firstEdgeId];
                secondEdge = &PM.Cells1D[secondEdgeId];
            }

            unsigned int newVertice = firstEdgeId == firstEdge->Vertices[0] ? firstEdge->Vertices[1] : firstEdge->Vertices[0];
            newCell2DVertices.push_back(newVertice);

            firstEdge->NearCells2D.push_back(newCell2D.Id);
            secondEdge->NearCells2D.push_back(newCell2D.Id);

            newCell2DEdges.push_back(firstEdgeId);
            newCell2DEdges.push_back(secondEdgeId);
        }
        else
        {
            cell->NearCells2D.push_back(newCell2D.Id);
            newCell2DEdges.push_back(i);
        }
    }

    newCell2D.NumberVertices = newCell2DVertices.size();
    newCell2D.Vertices = newCell2DVertices;
    newCell2D.NumberEdges = newCell2DEdges.size();
    newCell2D.Edges = newCell2DEdges;

    PM.NumberCell2D++;
    PM.Cells2D.push_back(newCell2D);
}


void CreateNewCells(PolygonalMesh &PM,
                    const DFNLibrary::Trace &t,
                    unsigned int &idCell0D,
                    unsigned int &idCell1D,
                    unsigned int &idCell2D,
                    const bool &pass,
                    const double &tol)
{
    Eigen::Vector3d p_r = t.EndpointsCoordinates.row(0);
    Eigen::Vector3d t_r = t.EndpointsCoordinates.row(1) - p_r.transpose();

    unsigned int indexCell2D = 0;

    std::list<unsigned int> cells2D;

    Eigen::Vector2d betaTemp = {};

    for(unsigned int i = 0; i < PM.Cells2D.size(); i++)
    {
        Cell2D *cell = &PM.Cells2D[i];

        if(!cell->IsOld && CellContainsTrace(PM, cell, p_r, t_r, betaTemp, tol))
        {
            cells2D.push_back(cell->Id);
            break;
        }
    }

    // Finchè ho celle vicine
    while(cells2D.size() != 0)
    {
        indexCell2D = cells2D.front();

        if(PM.Cells2D[indexCell2D].IsOld)
        {
            cells2D.erase(cells2D.begin());
            continue;
        }

        Cell2D *cell2D = &PM.Cells2D[indexCell2D];

        std::vector<unsigned int> indexIntersectedCell1D = {};

        Eigen::Vector2d beta = {};

        if(!pass && !CellContainsTrace(PM, cell2D, p_r, t_r, beta, tol))
        {
            cell2D->IsOld = true;
            cells2D.erase(cells2D.begin());
            CreateLastCell(PM, idCell2D, cell2D);
            continue;
        }

        std::vector<Cell1D> newEdges = {};

        std::vector<Cell0D> newVertices = {};

        unsigned int indexNewEdge = 0;

        bool first = true;

        std::vector<unsigned int> firstCell2DVertices = {};
        std::vector<unsigned int> secondCell2DVertices = {};

        std::vector<unsigned int> firstCell2DEdges = {};
        std::vector<unsigned int> secondCell2DEdges = {};

        bool repeat = false;

        for(unsigned int i = 0; i < cell2D->Edges.size(); i++)
        {
            if(repeat)
            {
                repeat = false;
                continue;
            }

            Cell1D *cell1D = &PM.Cells1D[cell2D->Edges[i]];

            Cell0D firstCell0D;
            Cell0D secondCell0D;

            // Salvo gli estremi della cella 1D in senso antiorario
            for(unsigned int j = 0; j < cell2D->Vertices.size(); j++)
            {
                int next = j == cell2D->Vertices.size() - 1 ? 0 : j + 1;

                int firstVertice = cell2D->Vertices[j];

                auto it = find(cell1D->Vertices.begin(), cell1D->Vertices.end(), firstVertice);

                if (it == cell1D->Vertices.end())
                    continue;

                int secondVertice = cell2D->Vertices[next];

                it = find(cell1D->Vertices.begin(), cell1D->Vertices.end(), secondVertice);

                if(it == cell1D->Vertices.end())
                    continue;

                firstCell0D = PM.Cells0D[firstVertice];
                secondCell0D = PM.Cells0D[secondVertice];
                break;
            }

            Eigen::Vector3d p_s = firstCell0D.Coordinates;
            Eigen::Vector3d t_s = secondCell0D.Coordinates - p_s;

            Eigen::Vector3d coordinates = {};

            if(first)
                firstCell2DVertices.push_back(firstCell0D.Id);
            else
                secondCell2DVertices.push_back(firstCell0D.Id);

            if(IntersectionCellTrace(p_s, t_s, p_r, t_r, coordinates, tol))
            {
                // Mantengo l'ordine delle celle intersecate
                auto it = find(indexIntersectedCell1D.begin(), indexIntersectedCell1D.end(), cell1D->Id);

                if(it == indexIntersectedCell1D.end())
                    indexIntersectedCell1D.push_back(cell1D->Id);

                unsigned int id;
                Cell0D newCell0D;
                Cell1D firstCell1D;
                Cell1D secondCell1D;

                if(AlreadyExists(coordinates, PM, id, tol))
                {
                    newCell0D = PM.Cells0D[id];
                    newVertices.push_back(newCell0D);

                    firstCell2DVertices.push_back(newCell0D.Id);
                    secondCell2DVertices.push_back(newCell0D.Id);

                    if(newCell0D.Id == secondCell0D.Id)
                    {
                        firstCell1D = PM.Cells1D[cell2D->Edges[i]];
                        secondCell1D = i == cell2D->Edges.size() - 1 ? PM.Cells1D[cell2D->Edges[0]] : PM.Cells1D[cell2D->Edges[i + 1]];

                        repeat = true;
                    }
                    else
                    {
                        unsigned int counter = 0;

                        Eigen::Vector2i firstCellVertices = {newCell0D.Id, firstCell0D.Id};
                        Eigen::Vector2i secondCellVertices = {secondCell0D.Id, newCell0D.Id};

                        for(Cell1D &cell : PM.Cells1D)
                        {
                            if(cell.Vertices == firstCellVertices)
                            {
                                firstCell1D = cell;
                                counter++;
                            }
                            else if(cell.Vertices == secondCellVertices)
                            {
                                secondCell1D = cell;
                                counter++;
                            }

                            if(counter == 2)
                                break;
                        }

                        cell1D->ReplacedBy = {firstCell1D.Id, secondCell1D.Id};
                    }
                }
                else
                {
                    newCell0D.Id = idCell0D++;
                    newCell0D.Coordinates = coordinates;

                    PM.NumberCell0D++;
                    PM.Cells0D.push_back(newCell0D);

                    newVertices.push_back(newCell0D);

                    firstCell2DVertices.push_back(newCell0D.Id);
                    secondCell2DVertices.push_back(newCell0D.Id);

                    cell1D->IsOld = true;

                    firstCell1D.Id = idCell1D++;
                    firstCell1D.Vertices = {firstCell0D.Id, newCell0D.Id};

                    newEdges.push_back(firstCell1D);

                    secondCell1D.Id = idCell1D++;
                    secondCell1D.Vertices = {newCell0D.Id, secondCell0D.Id};

                    cell1D->ReplacedBy = {firstCell1D.Id, secondCell1D.Id};

                    PM.NumberCell1D++;
                    PM.Cells1D.push_back(firstCell1D);

                    PM.NumberCell1D++;
                    PM.Cells1D.push_back(secondCell1D);

                    newEdges.push_back(secondCell1D);
                }

                if(first)
                {
                    firstCell2DEdges.push_back(firstCell1D.Id);

                    indexNewEdge = firstCell2DEdges.size(); // Salvo la posizione nella quale dovrò aggiungere la cella 1D della traccia
                    firstCell2DEdges.push_back(-1);

                    secondCell2DEdges.push_back(secondCell1D.Id);
                }
                else
                {
                    firstCell2DEdges.push_back(secondCell1D.Id);
                    secondCell2DEdges.push_back(firstCell1D.Id);
                }

                first = first ? false : true;
            }
            else
            {
                if(first)
                    firstCell2DEdges.push_back(cell1D->Id);
                else
                    secondCell2DEdges.push_back(cell1D->Id);
            }
        }

        cell2D->IsOld = true;

        CreateNewCellsByTrace(PM, idCell1D, idCell2D, newVertices, indexNewEdge,
                              firstCell2DVertices, secondCell2DVertices, firstCell2DEdges, secondCell2DEdges, tol);

        cells2D.erase(cells2D.begin());

        // Aggiungo alla lista le celle 2D se non le ho già intersecate
        for(unsigned int i : indexIntersectedCell1D)
        {
            for(unsigned int j : PM.Cells1D[i].NearCells2D)
            {
                if(!PM.Cells2D[j].IsOld)
                    cells2D.push_back(j);
            }
        }
    }
}

void GenerateMesh(PolygonalMesh &PM,
                  const DFNLibrary::Fracture &f,
                  const double &tol)
{

    unsigned int idCell2D = 0;

    unsigned int idCell0D = f.NumberVertices;
    unsigned int idCell1D = f.NumberVertices;

    // Memorizzo i dati della frattura come celle

    CreateFirstCell(PM, f, idCell2D);


    // Inizio a creare le celle nuove

    bool pass = true;

    for(const DFNLibrary::Trace &t : f.PassTraces)
    {
        if(t.IsOnEdge.at(f.Id))
            continue;

        CreateNewCells(PM, t, idCell0D, idCell1D, idCell2D, pass, tol);
    }

    pass = false;

    for(const DFNLibrary::Trace &t : f.NotPassTraces)
    {
        if(t.IsOnEdge.at(f.Id))
            continue;

        CreateNewCells(PM, t, idCell0D, idCell1D, idCell2D, pass, tol);
    }
}

}
