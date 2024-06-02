#include <iostream>
#include <fstream>
#include <sstream>
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
    // string fileName = "./DFN/"  + name + ".txt";

    string fileName = "C:/Users/user/OneDrive/Desktop/Progetto_PCS_2024/Project/src_test/TEST_data.txt";

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

    return 0;
}
