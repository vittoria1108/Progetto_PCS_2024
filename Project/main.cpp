#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "DFN.hpp"


using namespace std;
using namespace FractureLibrary;


int main(int argc, char ** argv)
{
    /*Vector3d normal1 = {1, 2, 3};
    Vector3d normal2 = {4, 5, 6};
    Vector3d t = normal1.cross(normal2);

    Matrix3d normalMatrix;
    normalMatrix.row(0) = normal1;
    normalMatrix.row(1) = normal2;
    normalMatrix.row(2) = t;
    normalMatrix = normalMatrix.transpose();

    for(unsigned int i = 0; i < 3; i++)
    {
        for(unsigned int j = 0; j < 3; j++)
        {
            cout << normalMatrix(i, j) << " ";
        }
        cout << endl;
    }*/

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
