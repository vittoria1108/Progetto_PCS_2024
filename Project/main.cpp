#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "DFN.hpp"


using namespace std;
using namespace FractureLibrary;


int main(int argc, char **argv)
{
    DFN dfn;
    if(argc == 1)
        return 2;
    istringstream str(argv[1]);
    string name;
    str >> name;
    string filename = "./DFN/" + name;

    if(!ImportFracture(filename, dfn))
    {
        return 1;
    }

    /*
    Vector4d aa = CalculatePlane(dfn.Fractures[1]);
    cout << "Coefficient of the plane" << endl;
    cout << "[ ";
    for (unsigned int i = 0; i<4; i++)
        cout << aa[i] << " ";
    cout << "] "<< endl;
    */
    for (unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        for (unsigned int j = 0; j< dfn.NumberFractures; j++)
        {
            if (i<j)
                CalculateTraces(dfn.Fractures[i], dfn.Fractures[j]);
        }
    }



    return 0;
}
