#include <iostream>
#include <Eigen/Eigen>
#include "DFN.hpp"
#include "PolygonalMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace FractureLibrary;

namespace PolygonalLibrary{

/*bool AlreadyExists(const Vector3d &coordinates, const PolygonalMesh &PM, unsigned int &id, const double &tol)
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
}*/

void ImportMesh(PolygonalMesh &PM, const Fracture &f, const double &tol)
{

    unsigned int idCell2D = 0;

    // Memorizzo i dati della frattura come celle

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

    unsigned int idCell0D = f.NumberVertices;
    unsigned int idCell1D = f.NumberVertices;

    // Inizio a creare le celle nuove

    for(const Trace &t : f.pTraces)
    {
        if(t.IsOnEdge.at(f.Id))
            continue;

        Vector3d p_r = t.EndpointsCoordinates.row(0);
        Vector3d t_r = t.EndpointsCoordinates.row(1) - p_r.transpose();

        unsigned int indexCell1D = 0;
        unsigned int indexCell2D = 0;

        vector<unsigned int> cells2D = PM.Cells1D[indexCell1D].NearCells2D;

        // QUANDO UNA TRACCIA INCONTRA DUE FRATTURE CONSIDERARE IL FATTO CHE UN VERTICE è IN COMUNE
        while(cells2D.size() != 0)
        {
            indexCell2D = cells2D[0];

            if(PM.Cells2D[indexCell2D].IsOld)
            {
                cells2D.erase(cells2D.begin());
                continue;
            }

            Vector2i vertices;
            unsigned int indexVertices = 0;
            unsigned int indexNewVertice = 0;

            Cell2D cell2D = PM.Cells2D[indexCell2D];

            bool first = true;

            vector<unsigned int> firstCell2DVertices = {};
            vector<unsigned int> secondCell2DVertices = {};

            vector<unsigned int> firstCell2DEdges = {};
            vector<unsigned int> secondCell2DEdges = {};

            for(unsigned int i : cell2D.Edges)
            {
                Cell1D *cell1D = &PM.Cells1D[i];

                Cell0D firstCell0D = PM.Cells0D[cell1D->Vertices[0]];
                Cell0D secondCell0D = PM.Cells0D[cell1D->Vertices[1]];

                Vector3d p_s = firstCell0D.Coordinates;
                Vector3d t_s = secondCell0D.Coordinates - p_s;

                Vector3d coordinates = {};

                if(first)
                    firstCell2DVertices.push_back(firstCell0D.Id);
                else
                    secondCell2DVertices.push_back(firstCell0D.Id);

                if(cell1D->IntersectsLine(p_s, t_s, p_r, t_r, coordinates, tol))
                {
                    Cell0D newCell0D;
                    newCell0D.Id = idCell0D++;
                    newCell0D.Coordinates = coordinates;

                    PM.NumberCell0D++;
                    PM.Cells0D.push_back(newCell0D);

                    vertices[indexVertices++] = newCell0D.Id;

                    PM.Cells1D[i].IsOld = true;

                    Cell1D newCell1;
                    newCell1.Id = idCell1D++;
                    newCell1.Vertices = {firstCell0D.Id, newCell0D.Id};

                    PM.NumberCell1D++;
                    PM.Cells1D.push_back(newCell1);

                    Cell1D newCell2;
                    newCell2.Id = idCell1D++;
                    newCell2.Vertices = {newCell0D.Id, secondCell0D.Id};

                    PM.NumberCell1D++;
                    PM.Cells1D.push_back(newCell2);

                    firstCell2DVertices.push_back(newCell0D.Id);
                    secondCell2DVertices.push_back(newCell0D.Id);

                    if(first)
                    {
                        firstCell2DEdges.push_back(newCell1.Id);
                        indexNewVertice = firstCell2DEdges.size(); // CONTROLLARE SE LA POSIZIONE è CORRETTA

                        secondCell2DEdges.push_back(newCell2.Id);
                    }
                    else
                    {
                        firstCell2DEdges.push_back(newCell2.Id);
                        secondCell2DEdges.push_back(newCell1.Id);
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

            Cell1D trace;
            trace.Id = idCell1D++;
            trace.Vertices = vertices;

            firstCell2DEdges.insert(firstCell2DEdges.begin() + indexNewVertice, trace.Id);
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
                PM.Cells1D[i].NearCells2D.push_back(firstCell2D.Id);
            }

            for(unsigned int i : secondCell2DEdges)
            {
                PM.Cells1D[i].NearCells2D.push_back(secondCell2D.Id);
            }

            trace.NearCells2D.push_back(firstCell2D.Id);
            trace.NearCells2D.push_back(secondCell2D.Id);

            PM.NumberCell1D++;
            PM.Cells1D.push_back(trace);

            cell2D.IsOld = true;
            cells2D.erase(cells2D.begin());
        }
    }
}

}
