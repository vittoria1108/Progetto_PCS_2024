#include <iostream>
#include <sstream>
#include <fstream>
#include "Utils.hpp"
#include <iomanip>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace DFNLibrary{

bool CompareTraces(const Trace &t1,
                   const Trace &t2)
{
    return t1.Length > t2.Length;
}


bool ImportFracture(const string &fileName,
                    DFN& dfn)
{

    ifstream file(fileName);

    if(file.fail())
        return false;


    string header;
    string line;

    getline(file, header);

    getline(file, line);

    istringstream nfractures(line);
    nfractures >> dfn.NumberFractures;

    dfn.Fractures.resize(dfn.NumberFractures);

    for(unsigned int i=0; i<dfn.NumberFractures; i++)
    {
        getline(file, header);

        getline(file, line, ';');
        istringstream idfractures(line);
        idfractures >> dfn.Fractures[i].Id;

        getline(file, line);
        istringstream nvertices(line);
        nvertices >> dfn.Fractures[i].NumberVertices;

        getline(file, header);

        dfn.Fractures[i].VerticesCoordinates.resize(3 , dfn.Fractures[i].NumberVertices);
        for(unsigned int j = 0; j < 3; j++)
        {
            getline(file, line);
            replace(line.begin(), line.end(), ';', ' ');
            istringstream converter(line);

            for(unsigned int k = 0; k < dfn.Fractures[i].NumberVertices; k++)
            {
                double coordinate;
                converter >> coordinate;
                dfn.Fractures[i].VerticesCoordinates(j,k) = coordinate;
            }
        }
        Vector3d barycentre = dfn.Fractures[i].VerticesCoordinates.rowwise().mean();
        dfn.Fractures[i].Barycentre = barycentre;
    }

    file.close();

    return true;
}


Vector4d CalculatePlane(const Fracture& fracture)
{
    Vector3d vector1 = fracture.VerticesCoordinates.col(1) - fracture.VerticesCoordinates.col(0);
    Vector3d vector2 = fracture.VerticesCoordinates.col(2) - fracture.VerticesCoordinates.col(0);

    Vector3d normal = vector1.cross(vector2);    // prodotto vettoriale

    Vector3d point = fracture.VerticesCoordinates.col(0);
    double d = normal.dot(point);     //normal[0]*point[0] + normal[1]*point[1] + normal[2]*point[2]

    Vector4d plane = {normal[0], normal[1], normal[2], d};

    return plane;   // restituisce i coefficienti dell'equazione del piano
}


double CalculateDistance(const Vector3d point1,
                         const Vector3d point2)
{
    double result = sqrt((point1-point2).transpose() * (point1-point2));   // distanza fra due punti in norma euclidea

    return result;
}


double CalculateR(const Fracture &f)
{
    double r = 0;
    for(unsigned int i = 0; i < f.NumberVertices; i++)
    {
        double newR = CalculateDistance(f.Barycentre, f.VerticesCoordinates.col(i));   // distanza fra il baricentro e tutti i vertici

        if(newR > r)
            r = newR;    // come raggio considero la distanza maggiore
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

    t_r = normal1.cross(normal2);  // direzione della retta di intersezione tr
    if(t_r.norm() < tol)
    {
        return false;  // se t=0 sono sullo stesso piano
    }


    // Se i piani non sono paralleli trovo la retta di intersezione
    Matrix3d normalMatrix;
    normalMatrix.row(0) = normal1;
    normalMatrix.row(1) = normal2;
    normalMatrix.row(2) = t_r;

    Vector3d constantTerms = {plane1[3], plane2[3], 0};

    p_r = normalMatrix.fullPivLu().solve(constantTerms);  // calcolo la fattorizzazione PALU per risolvere il sistema lineare

    return true;
}


bool IntersectionFractureLine(const Fracture &f,
                              const Vector3d &p_r,
                              const Vector3d &t_r,
                              Vector2d &beta,
                              const double &tol)
{
    // Controllo se la retta interseca le figure  (se le figure non sono tagliate dalla retta, non si intersecano)

    unsigned int counter = 0;

    for(unsigned int i = 0; i < f.NumberVertices; i++)
    {
        unsigned int next;     //per tenere conto anche dell'ultimo lato
        if (i == (f.NumberVertices -1))
            next = 0;
        else
            next = i+1;


        Vector3d p_s = f.VerticesCoordinates.col(i);
        Vector3d t_s = f.VerticesCoordinates.col(next) - f.VerticesCoordinates.col(i);  // direzione del segmento ts

        // segmento: p_s + aplha*t_s
        // retta: p_r + beta*t_r

        Vector3d prod = t_s.cross(t_r);


        if(prod.norm() > tol)   // controllo che t_s e t_r non siano parallele/coincidenti
        {
            double alpha = ((p_r-p_s).cross(t_r)).dot(prod)/(prod.dot(prod));
            // cout << fixed << setprecision(4) << "alpha1 " << alpha << endl;

            if(alpha > 0 && alpha < 1)   // combinazione convessa: se alpha sta tra 0 e 1 il segmento interseca effettivamente la retta t_r
            {
                double betaTemp = ((p_s-p_r).cross(t_s)).dot(-prod)/(prod.dot(prod));
                beta[counter] = betaTemp;
                counter++;

                if(counter==2)  // perché ho già trovato i due punti di intersezione
                {
                    sort(beta.begin(), beta.end());
                    return true;
                }
            }
        }
        else if((abs((p_s[0]-p_r[0]) / t_r[0] - (p_s[1] - p_r[1]) / t_r[1])) < tol && abs(((p_s[0]-p_r[0]) / t_r[0] - (p_s[2] - p_r[2] / t_r[2]))))
        {
            double betaTemp = ((p_s-p_r).dot(t_r))/(t_r.dot(t_r));
            beta[counter] = betaTemp;

            counter++;

            if(counter==2)  // perché ho già trovato i due punti di intersezione
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

    double barDistance = CalculateDistance(f1.Barycentre, f2.Barycentre);
    if(barDistance - (r1+r2) > tol)       // le due fratture non si intersecano
        return;


    // Controllo che i piani contenenti i due poligoni non siano paralleli
    Vector4d plane1 = CalculatePlane(f1);
    Vector4d plane2 = CalculatePlane(f2);
    Vector3d p_r;
    Vector3d t_r;

    if(!FindIntersectionLine(plane1, plane2, p_r, t_r, tol))
        return;

    // Controllo se la retta interseca le figure

    Vector2d beta_1 = {};

    if(!IntersectionFractureLine(f1, p_r, t_r, beta_1, tol))
        return;

    Vector2d beta_2 = {};

    if(!IntersectionFractureLine(f2, p_r, t_r, beta_2, tol))
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


    // passanti e non passanti
    if(beta_1[0] > beta_2[0] && beta_1[1] < beta_2[1])
    {
        f1.Tips[trace.Id] = false;      // passante = false
        f1.pTraces.push_back(trace);

        f2.Tips[trace.Id] = true;     // non passante = true
        f2.npTraces.push_back(trace);
    }
    else if(beta_2[0] > beta_1[0] && beta_2[1] < beta_1[1])
    {
        f1.Tips[trace.Id] = true;
        f1.npTraces.push_back(trace);

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
        f1.npTraces.push_back(trace);

        f2.Tips[trace.Id] = true;
        f2.npTraces.push_back(trace);
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
        tipsFile << f.Id << "; " << f.pTraces.size() + f.npTraces.size() << endl;

        tipsFile << "# TraceId; Tips; Length" << endl;

        for(const Trace &t : f.npTraces)
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


bool ReadDFN(const string &fileName,
             DFN &dfn,
             const double &tol)
{
    if(!ImportFracture(fileName,
                        dfn))
        return false;

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
        sort(dfn.Fractures[i].npTraces.begin(), dfn.Fractures[i].npTraces.end(), CompareTraces);
        sort(dfn.Fractures[i].pTraces.begin(), dfn.Fractures[i].pTraces.end(), CompareTraces);
    }

    return true;
}

}
