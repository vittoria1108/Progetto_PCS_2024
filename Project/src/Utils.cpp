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
    //cout << header << endl;       // # Number of Fractures

    getline(file, line);

    istringstream nFractures(line);
    nFractures >> dfn.NumberFractures;
    //cout << dfn.NumberFractures << endl;   // 3

    dfn.FracturesId.reserve(dfn.NumberFractures);

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        getline(file, header);
        //cout << header << endl;    // # FractureId; NumVertices

        getline(file, line, ';');
        istringstream idFractures(line);
        unsigned int id;
        idFractures >> id;
        dfn.FracturesId.push_back(id);
        //cout << "id: " << id << endl;  // id

        getline(file, line);
        istringstream nVertices(line);
        unsigned int nvert;
        nVertices >> nvert;
        dfn.NumberVertices.push_back(nvert);
        //cout << "numVertices: " << dfn.NumberVertices << endl;   // numVertices



        getline(file, header);
        //cout << header << endl;   // # Vertices

        vector<Vector3d> vertices;
        vertices.resize(nvert);
        for(unsigned int j = 0; j < 3; j++)
        {
            getline(file, line);
            replace(line.begin(), line.end(), ';', ' ');
            //cout << line << endl;
            istringstream converter(line);

            for(unsigned int k = 0; k < nvert; k++)
            {
                //double coordinate;
                converter >> vertices[k][j];
                //dfn.VerticesCoordinates[i][k][j] = coordinate;
            }

        }
        dfn.VerticesCoordinates[i] = vertices;



        /*
        for(unsigned int j = 0; j < dfn.NumberVertices; j++)
        {
            for(unsigned int k = 0; k < 3; k++)
            {
                cout << fixed << setprecision(16) << dfn.VerticesCoordinates[j][k] << endl;
            }
        }
        cout << endl;
        */
    }
    //cout << endl;

    /*
    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        cout << dfn.FracturesId[i] << " ";
    }
    cout << endl;
    */

    for(const auto& couple : dfn.VerticesCoordinates) //legge le coppie (chiave, vettori) nella mappa
    {
        cout<< couple.first << endl;
        for(const auto& vec : couple.second) //legge i vettori contenuti in couple.second
        {
            cout << "[ ";
            for(unsigned int c = 0; c < 3; c++)
                cout << fixed << setprecision(16) << vec(c) << " ";
            cout << "]" << endl;
        }
        /*
        METODO ALTERNATIVO per ottenere la stampa (o per ottenre separatamente id e vettore dei vertici separati)
        unsigned int id = couple.first;
        vector<Vector3d> vert = couple.second;
        cout << id << endl;
        for (unsigned int b = 0; b < vert.size(); b++)
        {
            cout << "[ ";
            for(unsigned int c = 0; c < 3; c++)
                cout << fixed << setprecision(16) << vert[b][c] << " ";
            cout << "]" << endl;

        }


        */
    }


    file.close();


    return true;



}

}
