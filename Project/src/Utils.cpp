#include <iostream>
#include <sstream>
#include <fstream>
#include "Utils.hpp"
#include <iomanip>


using namespace std;

namespace FractureLibrary{


bool ImportFracture(const string &filename,
                    DFN& dfn)
{

    ifstream file(filename);

    if(file.fail())
        return false;


    string header;
    string line;

    getline(file, header);
    cout << header << endl;       // # Number of Fractures

    getline(file, line);

    istringstream nfractures(line);
    nfractures >> dfn.NumberFractures;
    cout << dfn.NumberFractures << endl;   // 3


    for(unsigned int i=0; i<dfn.NumberFractures; i++)
    {
        getline(file, header);
        cout << header << endl;    // # FractureId; NumVertices

        getline(file, line, ';');
        istringstream idfractures(line);
        unsigned int id;
        idfractures >> id;
        dfn.FracturesId.push_back(id);
        cout << "id:" << " " << id << endl;  // id

        getline(file, line);
        istringstream nvertices(line);
        nvertices >> dfn.NumberVertices;
        cout << "numVertices:" << " " << dfn.NumberVertices << endl;   // numVertices

        getline(file, header);
        cout << header << endl;   // # Vertices


        dfn.VerticesCoordinates.reserve(dfn.NumberVertices);

        for(unsigned int j = 0; j < 3; j++)
        {
            getline(file, line);
            replace(line.begin(), line.end(), ';', ' ');
            istringstream converter(line);

            for(unsigned int k = 0; k < dfn.NumberVertices; k++)
            {
                double coordinate;
                converter >> coordinate;
                dfn.VerticesCoordinates[k][j] = coordinate;
            }
        }

        for(unsigned int k = 0; k < dfn.NumberVertices; k++)
        {
            for(unsigned int j = 0; j < 3; j++)
            {
                cout << scientific << setprecision(16) << dfn.VerticesCoordinates[k][j] << endl;
            }
        }
        cout << endl;
    }
    cout << endl;


    file.close();


    return true;

}

}
