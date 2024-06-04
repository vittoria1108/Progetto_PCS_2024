#include <iostream>
#include <Eigen/Eigen>
#include "DFN.hpp"
#include "PolygonalMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace DFNLibrary;

namespace PolygonalLibrary{

void ImportMesh(PolygonalMesh &PM,const Fracture &f)
{
    PM.NumberCell0D = f.NumberVertices;
    PM.NumberCell1D = f.NumberVertices;
    PM.NumberCell2D = 1;

    //resize dei vettori
    PM.IdCell0D.resize(PM.NumberCell0D);
    PM.CoordinatesCell0D.resize(PM.NumberCell0D);
    PM.IdCell1D.resize(PM.NumberCell1D);
    PM.Cell1DOld.resize(PM.NumberCell1D);
    PM.VerticesCell1D.resize(PM.NumberCell1D);
    PM.VerticesCell2D.resize(PM.NumberCell2D); //in realtà è 1 non so se ha senso
    PM.VerticesCell2D[0].resize(PM.NumberCell0D);
    PM.Cell2DOld.resize(PM.NumberCell2D); //in realtà è 1 non so se ha senso

    PM.Cell2DOld[0] = false; //per ora è tutto "nuovo"

    vector<unsigned int> vertices = {};
    for (unsigned int i=0; i< PM.NumberCell0D; i++)
    {
        unsigned int next = i+1;
        if (i == PM.NumberCell0D-1)
            next = 0;
        PM.IdCell0D[i] = i;
        PM.CoordinatesCell0D[i] = f.VerticesCoordinates.col(i);
        PM.IdCell1D[i] = i;
        PM.VerticesCell1D[i] = {i, next};
        PM.Cell1DOld[i] = false; //per ora è tutto "nuovo"
        PM.VerticesCell2D[0][i] = i;
    }

    cout << "Num0" << PM.NumberCell0D << endl;
    cout << "Num1" << PM.NumberCell1D << endl;
    cout << "Num2" << PM.NumberCell2D << endl;
}

};
