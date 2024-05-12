#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "DFN.hpp"


using namespace std;
using namespace FractureLibrary;


int main()
{
    DFN dfn;
    string filename = "DFN/FR3_data.txt";

    if(!ImportFracture(filename,
                        dfn))
    {
        return 1;
    }

    return 0;
}
