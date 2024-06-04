#include <iostream>
#include <fstream>
#include <sstream>
#include "src/Utils.hpp"
#include "src/DFN.hpp"
#include "src_mesh/PolygonalMesh.hpp"

using namespace std;
using namespace Eigen;
using namespace DFNLibrary;
using namespace PolygonalLibrary;


int main(int argc, char **argv)
{
    DFN dfn;
    if(argc == 1)
        return 1;

    istringstream str(argv[1]);
    string name;
    str >> name;
    double tol = 10 * numeric_limits<double>::epsilon(); //tolleranza di default

    if(argc == 3)//se viene inserita in command line la tolleranza
    {
        istringstream conv(argv[2]);
        double tol_input;
        conv >> tol_input;
        tol = max(tol_input, 10 * numeric_limits<double>::epsilon()); //scelgo la tol più alta tra default e input
    }

     string filename = "./DFN/" + name +".txt";
    if(!ReadDFN(filename,
                dfn,
                tol))
    {
        return 2;
    }

    string outputTracesFile = "Traces_" + name + ".txt";
    string outputTipsFile = "Tips_" + name + ".txt";

    WriteOutputFiles(outputTracesFile,
                     outputTipsFile,
                     dfn);
    /*
    PolygonalMesh PM;
    ImportMesh(PM , dfn.Fractures[1]);
    */

    return 0;
}
