#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "DFN.hpp"


using namespace std;
using namespace FractureLibrary;


int main(int argc, char ** argv)
{
    DFN dfn;

    if(argc == 1)
        return 1;

    istringstream str(argv[1]);
    string name;
    str >> name;
    string fileName = "./DFN/"  + name + ".txt";

    if(!ImportFracture(fileName,
                        dfn))
    {
        return 2;
    }

    string outputTracesFile = "Traces_" + name + ".txt";
    string outputTipsFile = "Tips_" + name + ".txt";

    WriteOutputFiles(outputTracesFile,
                     outputTipsFile,
                     dfn);

    return 0;
}
