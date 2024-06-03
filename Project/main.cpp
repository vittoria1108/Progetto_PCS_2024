#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "dfn.hpp"
#include <iomanip>


using namespace std;
using namespace DFNLibrary;


int main(int argc, char ** argv)
{
    if(argc == 1)
        return 1;

    istringstream str(argv[1]);
    string name;
    str >> name;
    string fileName = "./DFN/"  + name + ".txt";

    double tol = 10 * numeric_limits<double>::epsilon();   //tolleranza di default
    if(argc == 3)  //se viene inserita in command line la tolleranza
    {
        istringstream conv(argv[2]);
        double tol_input;
        conv >> tol_input;
        tol = max(tol_input, 10 * numeric_limits<double>::epsilon()); //scelgo la tol più alta tra default e input
    }

    DFN dfn;

    if(!ReadDFN(fileName,
                dfn,
                tol))
    {
        return 2;
    }
    else
    {
        cout << "Numero di fratture:" << " " << dfn.NumberFractures << endl;
        cout << endl;

        for(unsigned int i=0; i<dfn.NumberFractures; i++)
        {
            cout << "Id:" << " " << dfn.Fractures[i].Id << endl;
            cout << "Vertices:" << " " << dfn.Fractures[i].NumberVertices << endl;

            for(unsigned int j = 0; j < 3; j++)
            {
                cout << "[";
                for(unsigned int k = 0; k < dfn.Fractures[i].NumberVertices; k++)
                {
                    cout << scientific << setprecision(16) << dfn.Fractures[i].VerticesCoordinates(j,k) << " ";
                }
                cout << "]" << endl;
            }
            cout << endl;
        }
    }


    string outputTracesFile = "Traces_" + name + ".txt";
    string outputTipsFile = "Tips_" + name + ".txt";

    WriteOutputFiles(outputTracesFile,
                     outputTipsFile,
                     dfn);


    return 0;
}

