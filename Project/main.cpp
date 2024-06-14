#include <iostream>
#include <fstream>
#include <sstream>
#include "DFN.hpp"
#include "src_mesh/PolygonalMesh.hpp"


using namespace std;
using namespace FractureLibrary;
using namespace PolygonalLibrary;


int main(int argc, char ** argv)
{
    /*FractureLibrary::DFN dfn;

    if(argc == 1)
        return 1;

    std::istringstream str(argv[1]);
    std::string name;
    str >> name;
    std::string fileName = "./DFN/"  + name + ".txt";

    double tol = 10 * numeric_limits<double>::epsilon(); //tolleranza di default

    if(argc == 3)//se viene inserita in command line la tolleranza
    {
        std::istringstream conv(argv[2]);
        double tol_input;
        conv >> tol_input;
        tol = max(tol_input, 10 * numeric_limits<double>::epsilon()); //scelgo la tol più alta tra default e input
    }

    if(!ReadDFN(fileName,
                dfn,
                tol))
    {
        return 2;
    }

    std::string outputTracesFile = "Traces_" + name + ".txt";
    std::string outputTipsFile = "Tips_" + name + ".txt";

    WriteOutputFiles(outputTracesFile,
                     outputTipsFile,
                     dfn);

    std::vector<PolygonalLibrary::PolygonalMesh> allPM;

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        PolygonalLibrary::PolygonalMesh PM;
        ImportMesh(PM, dfn.Fractures[i], tol);

        allPM.push_back(PM);
    }*/

    Fracture f;
    f.Id = 0;
    f.NumberVertices = 4;

    MatrixXd vertices(3, 4);
    vertices << 0, 18, 18, 0,
                0, 0, 12, 12,
                0, 0, 0, 0;

    f.VerticesCoordinates = vertices;

    Trace t1;
    t1.Id = 0;
    t1.FracturesIds = {0, 1};
    t1.IsOnEdge[0] = false;

    MatrixXd coordinates1(2, 3);
    coordinates1 << 0, 6, 0,
                   18, 6, 0;

    t1.EndpointsCoordinates = coordinates1;

    Trace t2;
    t2.Id = 1;
    t2.FracturesIds = {0, 2};
    t2.IsOnEdge[0] = false;

    MatrixXd coordinates2(2, 3);
    coordinates2 << 9, 12, 0,
        9, 0, 0;

    t2.EndpointsCoordinates = coordinates2;

    Trace t3;
    t3.Id = 2;
    t3.FracturesIds = {0, 2};
    t3.IsOnEdge[0] = false;

    MatrixXd coordinates3(2, 3);
    coordinates3 << 10, 12, 0,
        10, 4, 0;

    t3.EndpointsCoordinates = coordinates3;

    Trace t4;
    t4.Id = 3;
    t4.FracturesIds = {0, 2};
    t4.IsOnEdge[0] = false;

    MatrixXd coordinates4(2, 3);
    coordinates4 << 12, 9, 0,
        7, 4, 0;

    t4.EndpointsCoordinates = coordinates4;

    Trace t5;
    t5.Id = 4;
    t5.FracturesIds = {0, 2};
    t5.IsOnEdge[0] = false;

    MatrixXd coordinates5(2, 3);
    coordinates5 << 0, 8, 0,
        5, 8, 0;

    t5.EndpointsCoordinates = coordinates5;

    f.pTraces.push_back(t1);
    f.pTraces.push_back(t2);

    f.nTraces.push_back(t3);
    f.nTraces.push_back(t4);
    f.nTraces.push_back(t5);


    PolygonalLibrary::PolygonalMesh PM;
    ImportMesh(PM, f, 1);

    return 0;
}
