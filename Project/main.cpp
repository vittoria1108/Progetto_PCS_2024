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
        return 2;

    istringstream str(argv[1]);
    string name;
    str >> name;
    string filename = "./DFN/"  + name;

    if(!ImportFracture(filename,
                        dfn))
    {
        return 1;
    }

    return 0;
}
