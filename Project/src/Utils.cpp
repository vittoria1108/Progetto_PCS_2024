#include <iostream>
#include <sstream>
#include <fstream>
#include "Utils.hpp"

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

    dfn.FracturesId.reserve(dfn.NumberFractures);

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

        getline(file, header);
        getline(file, header);
        getline(file, header);

        cout << endl;

    }
    cout << endl;

    for(unsigned int i=0; i<dfn.NumberFractures; i++)
    {
        cout << dfn.FracturesId[i] << " ";
    }
    cout << endl;





    file.close();


    return true;



}

}
