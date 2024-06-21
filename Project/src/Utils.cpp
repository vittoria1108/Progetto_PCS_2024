#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <DFN.hpp>
#include <Eigen/Eigen>


using namespace std;
using namespace Eigen;


namespace DFNLibrary{

double CalculateSquareDistance(const Eigen::Vector3d &point1,
                               const Eigen::Vector3d &point2)
{
    double result = (point1 - point2).transpose() * (point1 - point2);
    return result;
}

bool FindIntersectionLine(const Eigen::Vector4d &plane1,
                          const Eigen::Vector4d &plane2,
                          Eigen::Vector3d &p_r,
                          Eigen::Vector3d &t_r,
                          const double &tol)
{
    Eigen::Vector3d normal1 = {plane1[0], plane1[1], plane1[2]};
    Eigen::Vector3d normal2 = {plane2[0], plane2[1], plane2[2]};

    t_r = normal1.cross(normal2);

    if(t_r.norm() < tol) // I piani sono paralleli
    {
        return false;
    }

    // Se non sono paralleli trovo la retta di intersezione

    Eigen::Matrix3d normalMatrix;
    normalMatrix.row(0) = normal1;
    normalMatrix.row(1) = normal2;
    normalMatrix.row(2) = t_r;

    Eigen::Vector3d constantTerms = {plane1[3], plane2[3], 0};
    p_r = normalMatrix.fullPivLu().solve(constantTerms);

    return true;
}

bool CompareTraces(const Trace &t1,
                   const Trace &t2)
{
    return t1.Length > t2.Length;
}


bool ImportFracture(const std::string &fileName,
                    DFN &dfn)
{

    std::ifstream file(fileName);

    if(file.fail())
    {
        cerr << "Input file not found";
        return false;
    }

    std::string header;
    std::string line;

    getline(file, header);    // # Number of Fractures

    getline(file, line);

    std::istringstream nFractures(line);
    nFractures >> dfn.NumberFractures;
    dfn.Fractures.resize(dfn.NumberFractures);

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        getline(file, header);  // # FractureId; NumVertices

        getline(file, line, ';');
        std::istringstream idFractures(line);
        idFractures >> dfn.Fractures[i].Id;

        getline(file, line);
        std::istringstream nVertices(line);
        nVertices >> dfn.Fractures[i].NumberVertices;  // numVertices

        dfn.Fractures[i].VerticesCoordinates.resize(3, dfn.Fractures[i].NumberVertices);

        getline(file, header);  // # Vertices

        for(unsigned int j = 0; j < 3; j++)
        {
            getline(file, line);
            replace(line.begin(), line.end(), ';', ' ');
            std::istringstream converter(line);

            for(unsigned int k = 0; k < dfn.Fractures[i].NumberVertices; k++)
            {
                converter >> dfn.Fractures[i].VerticesCoordinates(j, k);
            }
        }

        Eigen::Vector3d barycentre = dfn.Fractures[i].VerticesCoordinates.rowwise().mean();
        dfn.Fractures[i].Barycentre = barycentre;
    }

    file.close();

    return true;
}

void CalculateTraces(DFN &dfn,
                     Fracture &f1,
                     Fracture &f2,
                     unsigned int &id,
                     const double &tol)
{
    // Controllo che le sfere che contengono i due poligoni non si intersecano

    double r1 = f1.CalculateR();
    double r2 = f2.CalculateR();

    double barDistances = CalculateSquareDistance(f1.Barycentre, f2.Barycentre);

    if(barDistances - ((r1 + r2) * (r1 + r2)) > tol)
        return;

    // Calcolo i piani contenenti i poligoni

    Eigen::Vector4d plane1 = f1.CalculatePlane();
    Eigen::Vector4d plane2 = f2.CalculatePlane();

    Eigen::Vector3d p_r;
    Eigen::Vector3d t_r;

    Eigen::Vector2d beta_1 = {};
    Eigen::Vector2d beta_2 = {};

    std::map<unsigned int, bool> isOnEdge;

    isOnEdge[f1.Id] = false;
    isOnEdge[f2.Id] = false;

    if(FindIntersectionLine(plane1, plane2, p_r, t_r, tol)) // Controllo se riesco a trovare la retta di intersezione tra i piani
    {                                                       // (se non sono paralleli)
        // Controllo se la retta interseca le figure

        if(!f1.IntersectsLine(p_r, t_r, beta_1, isOnEdge, tol))
            return;

        if(!f2.IntersectsLine(p_r, t_r, beta_2, isOnEdge, tol))
            return;

        if(beta_1[1] < beta_2[0])
            return;

        if(beta_2[1] < beta_1[0])
            return;
    }
    else    // Se i piani sono paralleli, controllo se effettivamente coincidono e, in questo caso,
    {       // se le fratture hanno dei lati che si intersecano
        if(f1.IsInPlane(plane2, tol) && f1.IntersectsEdges(f2, beta_1, beta_2, p_r, t_r, tol))
        {
            isOnEdge[f1.Id] = true;
            isOnEdge[f2.Id] = true;
        }
        else
            return;
    }


    // Dopo aver effettuato tutti i controlli so che la traccia esiste, quindi la creo

    Trace trace;
    trace.Id = id++;
    trace.FracturesIds = {f1.Id, f2.Id};

    Eigen::MatrixXd endPoints(2, 3);

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

    double squareLength = CalculateSquareDistance(endPoints.row(0), endPoints.row(1));

    if(abs(squareLength) < tol)
        return;

    trace.Length = sqrt(squareLength);

    trace.IsOnEdge = isOnEdge;

    if(abs(beta_1[0] - beta_2[0]) < tol && abs(beta_1[1] - beta_2[1]) < tol)
    {
        f1.Tips[trace.Id] = false;
        f1.PassTraces.push_back(trace);

        f2.Tips[trace.Id] = false;
        f2.PassTraces.push_back(trace);
    }
    else if(beta_1[0] > beta_2[0] && beta_1[1] < beta_2[1])
    {
        f1.Tips[trace.Id] = false;
        f1.PassTraces.push_back(trace);

        f2.Tips[trace.Id] = true;
        f2.NotPassTraces.push_back(trace);
    }
    else if(beta_2[0] > beta_1[0] && beta_2[1] < beta_1[1])
    {
        f1.Tips[trace.Id] = true;
        f1.NotPassTraces.push_back(trace);

        f2.Tips[trace.Id] = false;
        f2.PassTraces.push_back(trace);
    }
    else
    {
        f1.Tips[trace.Id] = true;
        f1.NotPassTraces.push_back(trace);

        f2.Tips[trace.Id] = true;
        f2.NotPassTraces.push_back(trace);
    }

    dfn.Traces.push_back(trace);
}

bool ReadDFN(const std::string &fileName,
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
        sort(dfn.Fractures[i].NotPassTraces.begin(), dfn.Fractures[i].NotPassTraces.end(), CompareTraces);
        sort(dfn.Fractures[i].PassTraces.begin(), dfn.Fractures[i].PassTraces.end(), CompareTraces);
    }

    return true;
}

void WriteOutputFiles(const std::string &outputTracesFile,
                      const std::string &outputTipsFile,
                      const DFN &dfn)
{
    std::ofstream tracesFile(outputTracesFile);

    tracesFile << "# Number of Traces" << std::endl;
    tracesFile << dfn.NumberTraces << std::endl;

    tracesFile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << std::endl;

    for(Trace trace : dfn.Traces)
    {
        tracesFile << trace.Id << "; " << trace.FracturesIds[0] << "; " << trace.FracturesIds[1] << "; ";

        for(unsigned int i = 0; i < 2; i++)
        {
            for(unsigned int j = 0; j < 3; j++)
            {
                if(i != 1 || j != 2)
                    tracesFile << setprecision(16) << scientific << trace.EndpointsCoordinates(i, j) << "; ";
            }
        }

        tracesFile << setprecision(16) << scientific << trace.EndpointsCoordinates(1, 2) << std::endl;
        tracesFile << std::endl;
    }

    tracesFile.close();

    std::ofstream tipsFile(outputTipsFile);

    for(const Fracture &f : dfn.Fractures)
    {
        unsigned int totalSize = f.PassTraces.size() + f.NotPassTraces.size();

        tipsFile << "# FractureId; NumTraces" << std::endl;
        tipsFile << f.Id << "; " << totalSize << std::endl;

        if(totalSize != 0)
            tipsFile << "# TraceId; Tips; Length" << std::endl;

        for(const Trace &t : f.PassTraces)
        {
            tipsFile << t.Id << "; false; " << setprecision(16) << scientific << t.Length << std::endl;
        }

        for(const Trace &t : f.NotPassTraces)
        {
            tipsFile << t.Id << "; true; " << setprecision(16) << scientific << t.Length << std::endl;
        }

        tipsFile << std::endl;
    }

    tipsFile.close();
}

}
