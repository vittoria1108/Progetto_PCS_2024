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
    DFN dfn;

    if(argc == 1)
        return 1;

    istringstream str(argv[1]);
    string name;
    str >> name;
    string fileName = "./DFN/"  + name + ".txt";

    double tol = 10 * numeric_limits<double>::epsilon(); //tolleranza di default

    if(argc == 3)//se viene inserita in command line la tolleranza
    {
        istringstream conv(argv[2]);
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

    string outputTracesFile = "Traces_" + name + ".txt";
    string outputTipsFile = "Tips_" + name + ".txt";

    WriteOutputFiles(outputTracesFile,
                     outputTipsFile,
                     dfn);

    vector<PolygonalMesh> allPM;

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        PolygonalMesh PM;
        ImportMesh(PM, dfn.Fractures[i], tol);

        allPM.push_back(PM);
    }

    return 0;
}
