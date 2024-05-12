#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "Utils.hpp"
#include "Eigen/Eigen"


using namespace std;
using namespace Eigen;


namespace FractureLibrary{

bool ImportFracture(const string &fileName,
                    DFN &dfn)
{

    ifstream file(fileName);

    if(file.fail())
        return false;

    string header;
    string line;

    getline(file, header);
    cout << header << endl;       // # Number of Fractures

    getline(file, line);

    istringstream nFractures(line);
    nFractures >> dfn.NumberFractures;
    cout << dfn.NumberFractures << endl;   // 3

    dfn.FracturesId.reserve(dfn.NumberFractures);

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        getline(file, header);
        cout << header << endl;    // # FractureId; NumVertices

        getline(file, line, ';');
        istringstream idFractures(line);
        unsigned int id;
        idFractures >> id;
        dfn.FracturesId.push_back(id);
        cout << "id: " << id << endl;  // id

        getline(file, line);
        istringstream nVertices(line);
        nVertices >> dfn.NumberVertices;
        cout << "numVertices: " << dfn.NumberVertices << endl;   // numVertices

        getline(file, header);
        cout << header << endl;   // # Vertices

        vector<Vector3d> vertices;
        vertices.resize(dfn.NumberVertices);

        for(unsigned int j = 0; j < 3; j++)
        {
            getline(file, line);
            replace(line.begin(), line.end(), ';', ' ');
            cout << line << endl;
            istringstream converter(line);

            for(unsigned int k = 0; k < dfn.NumberVertices; k++)
            {
                double coordinate;
                converter >> coordinate;
                vertices[k][j] = coordinate;
            }
        }

        for(unsigned int k = 0; k < dfn.NumberVertices; k++)
        {
            dfn.FracturesCoordinates.push_back(vertices[k]);
        }

        for(unsigned int j = 0; j < dfn.NumberVertices; j++)
        {
            for(unsigned int k = 0; k < 3; k++)
            {
                cout << fixed << setprecision(16) << dfn.FracturesCoordinates[j][k] << endl;
            }
        }
        cout << endl;

    }
    cout << endl;

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        cout << dfn.FracturesId[i] << " ";
    }
    cout << endl;





    file.close();


    return true;



}

}
