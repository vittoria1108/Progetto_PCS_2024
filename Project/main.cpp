#include <iostream>
#include <fstream>
#include <sstream>
#include <DFN.hpp>
#include "src_mesh/PolygonalMesh.hpp"


using namespace std;
using namespace FractureLibrary;
using namespace PolygonalLibrary;


int main(int argc, char ** argv)
{
    FractureLibrary::DFN dfn;

    if(argc == 1) // Non ho passato alcun file
        return 1;

    std::istringstream str(argv[1]);
    std::string name;
    str >> name;
    std::string fileName = "./DFN/"  + name + ".txt";

    double tol = 10 * numeric_limits<double>::epsilon(); // Tolleranza di default

    if(argc == 3) // Se viene inserita in command line la tolleranza
    {
        std::istringstream conv(argv[2]);
        double tol_input;
        conv >> tol_input;
        tol = max(tol_input, 10 * numeric_limits<double>::epsilon()); // Scelgo la tol più alta tra default e input
    }

    if(!ReadDFN(fileName,
                dfn,
                tol))
    {
        return 2;
    }

    std::string outputTracesFile = "Traces_" + name + ".txt";
    std::string outputTipsFile = "Tips_" + name + ".txt";

    WriteOutputFiles(outputTracesFile,
                     outputTipsFile,
                     dfn);


    /*unsigned int num_iter = 10;

    double average = 0;

    for(unsigned int k = 0; k < num_iter; k++)
    {
        std::chrono::steady_clock::time_point t_begin = std::chrono::steady_clock::now();
        ReadDFN(fileName,
                dfn,
                tol);
        std::chrono::steady_clock::time_point t_end = std::chrono::steady_clock::now();

        double duration = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_begin).count();

        average += duration;
        std::cout << duration << std::endl;
    }

    std::cout << std::endl;
    average = average/num_iter;

    std::cout << "Average time: " << average << std::endl;*/


    std::vector<PolygonalLibrary::PolygonalMesh> allPM;

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        PolygonalLibrary::PolygonalMesh PM;
        GenerateMesh(PM, dfn.Fractures[i], tol);

        allPM.push_back(PM);
    }

    unsigned int counter = 1;

    for(Cell2D &cell : allPM[0].Cells2D)
    {
        if(!cell.IsOld)
        {
            cout << "Cella " << counter++ << ":" << endl;

            for(unsigned int i : cell.Vertices)
            {
                cout << allPM[0].Cells0D[i].Coordinates[0] << "; " << allPM[0].Cells0D[i].Coordinates[1] << "; " << allPM[0].Cells0D[i].Coordinates[2] << endl;
            }
        }
    }

    return 0;
}
