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
    string filename = "C:/Users/user/OneDrive/Desktop/Progetto_PCS_2024/Project/DFN/FR3_data.txt";

    if(!ImportFracture(filename,
                        dfn))
    {
        return 1;
    }

    return 0;
}
