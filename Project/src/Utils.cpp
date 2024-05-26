#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "Utils.hpp"
#include "Eigen/Eigen"


using namespace std;
using namespace Eigen;


namespace FractureLibrary{

bool CompareTraces(const Trace &t1,
                   const Trace &t2)
{
    return t1.Length > t2.Length;
}

bool ImportFracture(const string &fileName,
                    DFN &dfn,
                    const double &tol)
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

    unsigned int traceId = 0;

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        for(unsigned int j = 0; j < dfn.NumberFractures; j++)
        {
            if(i < j)
            {
                CalculateTraces(dfn, dfn.Fractures[i], dfn.Fractures[j], traceId, tol);
            }
        }
    }

    dfn.NumberTraces = traceId;

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        sort(dfn.Fractures[i].nTraces.begin(), dfn.Fractures[i].nTraces.end(), CompareTraces);
        sort(dfn.Fractures[i].pTraces.begin(), dfn.Fractures[i].pTraces.end(), CompareTraces);
    }

    file.close();

    return true;

}

Vector4d CalculatePlane(const Fracture &f)
{
    Vector3d vector1 = f.VerticesCoordinates.col(1) - f.VerticesCoordinates.col(0);
    Vector3d vector2 = f.VerticesCoordinates.col(2) - f.VerticesCoordinates.col(0);

    Vector3d normal = vector1.cross(vector2);

    Vector3d point = f.VerticesCoordinates.col(0);
    double d = normal.dot(point);

    Vector4d plane = {normal[0], normal[1], normal[2], d};

    return plane;
}

double CalculateDistance(const Vector3d point1,
                         const Vector3d point2)
{
    double result = sqrt((point1 - point2).transpose() * (point1 - point2));
    return result;
}

double CalculateR(const Fracture &f)
{
    double r = 0;
    for(unsigned int i = 0; i < f.NumberVertices; i++)
    {
        double newR = CalculateDistance(f.Barycentre, f.VerticesCoordinates.col(i));

        if(newR > r)
            r = newR;
    }

    return r;
}

bool FindIntersectionLine(const Vector4d &plane1,
                          const Vector4d &plane2,
                          Vector3d &p_r,
                          Vector3d &t_r,
                          const double &tol)
{

    Vector3d normal1 = {plane1[0], plane1[1], plane1[2]};
    Vector3d normal2 = {plane2[0], plane2[1], plane2[2]};

    t_r = normal1.cross(normal2);
    if(t_r.norm() < tol)
    {
        return false;
    }

    // Se non sono paralleli trovo la retta di intersezione
    Matrix3d normalMatrix;
    normalMatrix.row(0) = normal1;
    normalMatrix.row(1) = normal2;
    normalMatrix.row(2) = t_r;

    Vector3d constantTerms = {plane1[3], plane2[3], 0};

    p_r = normalMatrix.fullPivLu().solve(constantTerms);

    return true;
}

bool IntersectionFractureLine(const Fracture &f,
                              const Vector3d &p_r,
                              const Vector3d &t_r,
                              Vector2d &beta,
                              const double &tol)
{
    unsigned int counter = 0;

    for(unsigned int i = 0; i < f.NumberVertices; i++)
    {
        unsigned int next; //per tenere conto anche dell'ultimo lato

        if (i == (f.NumberVertices -1))
            next = 0;
        else
            next = i+1;

        /* segmento s: Ps + alfa*ts
           retta    r: Pr + beta*tr */

        Vector3d p_s = f.VerticesCoordinates.col(i);
        Vector3d t_s = f.VerticesCoordinates.col(next) - f.VerticesCoordinates.col(i);
        Vector3d prod = t_s.cross(t_r);

        if (prod.norm() > tol) //controllo che non siano parallele
        {
            double alpha = (((p_r - p_s).cross(t_r)).dot(prod)) / (prod.dot(prod));

            if(alpha > 0 && alpha < 1)
            {
                double betaTemp = (((p_s - p_r).cross(t_s)).dot(-prod)) / (prod.dot(prod));
                beta[counter++] = betaTemp;

                if(counter == 2)
                {
                    sort(beta.begin(), beta.end());
                    return true;
                }

            }
        }
        // chiedere alla teora, c'è sicuramente una condizione più semplice
        else if ((abs((p_s[0]-p_r[0])/t_r[0] - (p_s[1]-p_r[1])/t_r[1]) < tol) && (abs((p_s[0]-p_r[0])/t_r[0] - (p_s[2]-p_r[2])/t_r[2])))//se ps appartiene alla retta t_r, quindi le due rette sono coincidenti e il lato giace sulla retta t_r
        {
            double betaTemp = ((p_s - p_r).dot(t_r))/(t_r).dot(t_r);
            beta[counter++] = betaTemp;

            if(counter == 2)
            {
                sort(beta.begin(), beta.end());
                return true;
            }

        }


    }


    return false;

}

void CalculateTraces(DFN &dfn,
                     Fracture &f1,
                     Fracture &f2,
                     unsigned int &id,
                     const double &tol)
{
    // Controllo che le circoferenze che contengono i due poligoni non si intersecano

    double r1 = CalculateR(f1);
    double r2 = CalculateR(f2);

    double barDistances = CalculateDistance(f1.Barycentre, f2.Barycentre);

    if(barDistances - (r1 + r2) > tol)
        return;

    // Controllo che i piani contenenti i due poligoni non siano paralleli
    Vector4d plane1 = CalculatePlane(f1);
    Vector4d plane2 = CalculatePlane(f2);
    Vector3d p_r;
    Vector3d t_r;

    if(!FindIntersectionLine(plane1, plane2, p_r, t_r,tol))
        return;

    // Controllo se la retta interseca le figure

    Vector2d beta_1 = {};

    if(!IntersectionFractureLine(f1, p_r, t_r, beta_1,tol))
        return;

    Vector2d beta_2 = {};

    if(!IntersectionFractureLine(f2, p_r, t_r, beta_2,tol))
        return;

    if(beta_1[1] < beta_2[0])
        return;

    if(beta_2[1] < beta_1[0])
        return;

    // Dopo aver effettuato tutti i controlli so che la traccia esiste, quindi la creo

    Trace trace;
    trace.Id = id++;
    trace.FracturesIds = {f1.Id, f2.Id};

    MatrixXd endPoints(2, 3);

    if(beta_1[0] < beta_2[0])
    {
        endPoints.row(0) = p_r + (beta_2[0] * t_r);
    }
    else
    {
        endPoints.row(0) = p_r + (beta_1[0] * t_r);
    }

    if(beta_1[1] > beta_2[1])
    {
        endPoints.row(1) = p_r + (beta_2[1] * t_r);
    }
    else
    {
        endPoints.row(1) = p_r + (beta_1[1] * t_r);
    }

    trace.EndpointsCoordinates = endPoints;

    double length = CalculateDistance(endPoints.row(0), endPoints.row(1));
    trace.Length = length;

    if(beta_1[0] > beta_2[0] && beta_1[1] < beta_2[1])
    {
        f1.Tips[trace.Id] = false;
        f1.pTraces.push_back(trace);

        f2.Tips[trace.Id] = true;
        f2.nTraces.push_back(trace);
    }
    else if(beta_2[0] > beta_1[0] && beta_2[1] < beta_1[1])
    {
        f1.Tips[trace.Id] = true;
        f1.nTraces.push_back(trace);

        f2.Tips[trace.Id] = false;
        f2.pTraces.push_back(trace);
    }
    else if(beta_1[0] == beta_2[0] && beta_1[1] == beta_2[1])
    {
        f1.Tips[trace.Id] = false;
        f1.pTraces.push_back(trace);

        f2.Tips[trace.Id] = false;
        f2.pTraces.push_back(trace);
    }
    else
    {
        f1.Tips[trace.Id] = true;
        f1.nTraces.push_back(trace);

        f2.Tips[trace.Id] = true;
        f2.nTraces.push_back(trace);
    }

    dfn.Traces.push_back(trace);
}

void WriteOutputFiles(const string &outputTracesFile,
                      const string &outputTipsFile,
                      const DFN &dfn)
{
    ofstream tracesFile(outputTracesFile);

    tracesFile << "# Number of Traces" << endl;
    tracesFile << dfn.NumberTraces << endl;

    tracesFile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

    for(Trace trace : dfn.Traces)
    {
        tracesFile << trace.Id << "; " << trace.FracturesIds[0] << "; " << trace.FracturesIds[1] << "; ";

        for(unsigned int i = 0; i < 2; i++)
        {
            for(unsigned int j = 0; j < 3; j++)
            {
                if(i != 1 || j != 2)
                    tracesFile << trace.EndpointsCoordinates(i, j) << "; ";
            }
        }

        tracesFile << trace.EndpointsCoordinates(1, 2) << endl;
    }

    tracesFile.close();

    ofstream tipsFile(outputTipsFile);

    for(const Fracture &f : dfn.Fractures)
    {
        tipsFile << "# FractureId; NumTraces" << endl;
        tipsFile << f.Id << "; " << f.pTraces.size() + f.nTraces.size() << endl;

        tipsFile << "# TraceId; Tips; Length" << endl;

        for(const Trace &t : f.nTraces)
        {
            tipsFile << t.Id << "; true; " << t.Length << endl;
        }

        for(const Trace &t : f.pTraces)
        {
            tipsFile << t.Id << "; false; " << t.Length << endl;
        }
        tipsFile << endl;
    }

    tipsFile.close();
}

}
