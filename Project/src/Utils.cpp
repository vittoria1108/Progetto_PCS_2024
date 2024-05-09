#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "Utils.hpp"


using namespace std;
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

        vector<double> xCoordinates = {};
        vector<double> yCoordinates = {};
        vector<double> zCoordinates = {};

        for(unsigned int j = 0; j < dfn.NumberVertices - 1; j++)
        {
            getline(file, line, ';');
            istringstream coordinate(line);
            double x;
            coordinate >> x;
            xCoordinates.push_back(x);
        }

        getline(file, line);
        istringstream coordinateX(line);
        double x;
        coordinateX >> x;
        xCoordinates.push_back(x);

        for(unsigned int j = 0; j < dfn.NumberVertices - 1; j++)
        {
            getline(file, line, ';');
            istringstream coordinate(line);
            double y;
            coordinate >> y;
            yCoordinates.push_back(y);
        }

        getline(file, line);
        istringstream coordinateY(line);
        double y;
        coordinateY >> y;
        yCoordinates.push_back(y);

        for(unsigned int j = 0; j < dfn.NumberVertices - 1; j++)
        {
            getline(file, line, ';');
            istringstream coordinate(line);
            double z;
            coordinate >> z;
            zCoordinates.push_back(z);
        }

        getline(file, line);
        istringstream coordinateZ(line);
        double z;
        coordinateZ >> z;
        zCoordinates.push_back(z);

        for(unsigned int j = 0; j < dfn.NumberVertices; j++)
        {
            Vector3d coordinates = {xCoordinates[j], yCoordinates[j], zCoordinates[j]};
            dfn.FracturesCoordinates.push_back(coordinates);
        }

        for(unsigned int j = 0; j < dfn.NumberVertices; j++)
        {
            for(unsigned int k = 0; k < 3; k++)
            {
                cout << setprecision(16) << dfn.FracturesCoordinates[j][k] << endl;
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
