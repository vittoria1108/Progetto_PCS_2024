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

    getline(file, header);    // # Number of Fractures

    getline(file, line);

    istringstream nFractures(line);
    nFractures >> dfn.NumberFractures;
    dfn.Fractures.resize(dfn.NumberFractures);

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        getline(file, header);  // # FractureId; NumVertices

        getline(file, line, ';');
        istringstream idFractures(line);
        idFractures >> dfn.Fractures[i].Id;

        getline(file, line);
        istringstream nVertices(line);
        nVertices >> dfn.Fractures[i].NumberVertices;  // numVertices

        dfn.Fractures[i].VerticesCoordinates.resize(3, dfn.Fractures[i].NumberVertices);

        getline(file, header);  // # Vertices

        for(unsigned int j = 0; j < 3; j++)
        {
            getline(file, line);
            replace(line.begin(), line.end(), ';', ' ');
            istringstream converter(line);

            for(unsigned int k = 0; k < dfn.Fractures[i].NumberVertices; k++)
            {
                converter >> dfn.Fractures[i].VerticesCoordinates(j, k);
            }
        }

        Vector3d barycentre = dfn.Fractures[i].VerticesCoordinates.rowwise().mean();
        dfn.Fractures[i].Barycentre = barycentre;
    }

    file.close();

    return true;

}

Vector4d CalculatePlane(const Fracture &f)
{
    Vector3d vector1 = f.VerticesCoordinates.col(1) - f.VerticesCoordinates.col(0);
    Vector3d vector2 = f.VerticesCoordinates.col(2) - f.VerticesCoordinates.col(0);

    Vector3d normal = vector1.cross(vector2);

    double norm = sqrt(normal.transpose() * normal);
    normal /= norm;

    Vector3d point = f.VerticesCoordinates.col(0);
    double d = -(normal[0]*point[0] + normal[1]*point[1] + normal[2]*point[2]);

    Vector4d plane = {normal[0], normal[1], normal[2], d};

    return plane;
}

double CalculateDistance(const Vector3d point1,
                         const Vector3d point2)
{
    double result = sqrt((point1 - point2).transpose() * (point1 - point2));
    return result;
}

bool CalculateTraces(const Fracture &f1,
                     const Fracture &f2)
{
    double r1 = 0;
    for(unsigned int i = 0; i < f1.NumberVertices; i++)
    {
        double newR = CalculateDistance(f1.Barycentre, f1.VerticesCoordinates.col(i));

        if(newR > r1)
            r1 = newR;
    }

    double r2 = 0;
    for(unsigned int i = 0; i < f2.NumberVertices; i++)
    {
        double newR = CalculateDistance(f2.Barycentre, f2.VerticesCoordinates.col(i));

        if(newR > r2)
            r2 = newR;
    }

    double barDistances = CalculateDistance(f1.Barycentre, f2.Barycentre);

    if(abs(barDistances - r1 - r2) > 10 * numeric_limits<double>::epsilon())
        return false;

    Vector4d plane1 = CalculatePlane(f1);
    Vector4d plane2 = CalculatePlane(f2);

    Vector3d normal1 = {plane1[0], plane1[1], plane1[2]};
    Vector3d normal2 = {plane2[0], plane2[1], plane2[2]};

    Vector3d t = normal1.cross(normal2);

    for(unsigned int i = 0; i < 3; i++)
    {
        if(abs(t[i]) > 10 * numeric_limits<double>::epsilon())
            break;

        return false;
    }

}

}
