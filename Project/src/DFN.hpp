#ifndef DFN_HPP
#define DFN_HPP

#include <iostream>
#include <iomanip>
#include <Eigen/Eigen>


using namespace std;
using namespace Eigen;


namespace FractureLibrary{

double CalculateSquareDistance(const Eigen::Vector3d point1,
                               const Eigen::Vector3d point2);

struct Trace{

    unsigned int Id;
    Eigen::Vector2i FracturesIds = {};
    Eigen::MatrixXd EndpointsCoordinates = {};
    double Length;

    std::map<unsigned int, bool> IsOnEdge = {};

};

struct Fracture{

    unsigned int Id = 0;
    unsigned int NumberVertices = 0;
    Eigen::MatrixXd VerticesCoordinates = {};
    Eigen::Vector3d Barycentre = {};
    std::map<unsigned int, bool> Tips = {};

    std::vector<Trace> nTraces = {};
    std::vector<Trace> pTraces = {};

    Eigen::Vector4d CalculatePlane()
    {
        Eigen::Vector3d vector1 = VerticesCoordinates.col(1) - VerticesCoordinates.col(0);
        Eigen::Vector3d vector2 = VerticesCoordinates.col(2) - VerticesCoordinates.col(0);

        Eigen::Vector3d normal = vector1.cross(vector2);

        Eigen::Vector3d point = VerticesCoordinates.col(0);
        double d = normal.dot(point);

        Eigen::Vector4d plane = {normal[0], normal[1], normal[2], d};

        return plane;
    }

    bool IsInPlane(const Eigen::Vector4d &plane,
                    const double &tol)
    {
        Eigen::Vector3d point = VerticesCoordinates.col(0);
        Eigen::Vector3d normal = {plane[0], plane[1], plane[2]};

        if(abs(point.dot(normal) - plane[3]) < tol)
            return true;

        return false;
    }

    double CalculateR()
    {
        double r = 0;
        for(unsigned int i = 0; i < NumberVertices; i++)
        {
            double newR = CalculateSquareDistance(Barycentre, VerticesCoordinates.col(i));

            if(newR > r)
                r = newR;
        }

        return sqrt(r);
    }

    bool IntersectsLine(const Eigen::Vector3d &p_r,
                        const Eigen::Vector3d &t_r,
                        Eigen::Vector2d &beta,
                        std::map<unsigned int, bool> &isOnEdge,
                        const double &tol)
    {
        unsigned int counter = 0;

        for(unsigned int i = 0; i < NumberVertices; i++)
        {
            unsigned int next; //per tenere conto anche dell'ultimo lato

            if (i == (NumberVertices - 1))
                next = 0;
            else
                next = i + 1;

            /* segmento s: Ps + alfa*ts
           retta    r: Pr + beta*tr */

            Eigen::Vector3d p_s = VerticesCoordinates.col(i);
            Eigen::Vector3d t_s = VerticesCoordinates.col(next) - VerticesCoordinates.col(i);
            Eigen::Vector3d prod = t_s.cross(t_r);

            if (prod.norm() > tol)
            {
                double alpha = (((p_r - p_s).cross(t_r)).dot(prod)) / (prod.dot(prod));

                if(alpha >= 0 && alpha < 1)
                {
                    double betaTemp = (((p_s - p_r).cross(t_s)).dot(-prod)) / (prod.dot(prod));
                    beta[counter++] = betaTemp;

                    if(counter == 2)
                    {
                        std::sort(beta.begin(), beta.end());
                        return true;
                    }
                }
            }
            else
            {
                bool sameLine = true;
                std::vector<unsigned int> indexBeta;
                double betaTemp = 0;
                double oldBeta = 0;

                // Se t_r[j] != 0 salvo l'indice, altrimenti controllo che p_s[j] == p_r[j]
                for(unsigned int j = 0; j < 3; j++)
                {
                    if(abs(t_r[j]) > tol)
                    {
                        indexBeta.push_back(j);
                    }
                    else if(abs(p_s[j] - p_r[j]) > tol)
                    {
                        sameLine = false;
                        break;
                    }
                }

                /* Se le coordinate dei punti sono uguali quando non riesco a calcolare le beta,
                    controllo che le beta calcolabili siano uguali*/
                if(sameLine)
                {
                    unsigned int index = indexBeta[0];
                    oldBeta = (p_s[index] - p_r[index]) / t_r[index];

                    indexBeta.erase(indexBeta.begin());

                    for(unsigned int j : indexBeta)
                    {
                        betaTemp = (p_s[j] - p_r[j]) / t_r[j];

                        if(abs(betaTemp - oldBeta) > tol)
                        {
                            sameLine = false;
                            break;
                        }

                        oldBeta = betaTemp;
                    }
                }

                if(sameLine)
                {
                    isOnEdge[Id] = true;

                    beta[counter++] = betaTemp;

                    if(counter == 2)
                    {
                        std::sort(beta.begin(), beta.end());
                        return true;
                    }
                }
            }
        }

        return false;
    }

    bool IntersectsEdges(Fracture &f,
                        Eigen::Vector2d &beta_1,
                        Eigen::Vector2d &beta_2,
                        Eigen::Vector3d &p_r,
                        Eigen::Vector3d &t_r,
                        const double &tol)
    {

        for(unsigned int i = 0; i < NumberVertices; i++)
        {
            unsigned int index = 0;

            unsigned int next1 = i == NumberVertices - 1 ? 0 : i + 1;

            // AB -> Lato prima frattura
            Eigen::Vector3d A = VerticesCoordinates.col(i);
            Eigen::Vector3d B = VerticesCoordinates.col(next1);

            Eigen::Vector3d line1 = B - A;

            p_r = A;
            t_r = line1;

            for(unsigned int k = 0; k < 3; k++)
            {
                if(abs(t_r[k]) > tol)
                {
                    index = k;
                    break;
                }
            }

            beta_1[0] = (A[index] - p_r[index]) / t_r[index];
            beta_1[1] = (B[index] - p_r[index]) / t_r[index];

            for(unsigned int j = 0; j < f.NumberVertices; j++)
            {
                bool sameLine = true;

                unsigned int next2 = j == f.NumberVertices - 1 ? 0 : j + 1;

                // CD -> Lato seconda frattura
                Eigen::Vector3d C = f.VerticesCoordinates.col(j);
                Eigen::Vector3d D = f.VerticesCoordinates.col(next2);

                beta_2[0] = (C[index] - p_r[index]) / t_r[index];
                beta_2[1] = (D[index] - p_r[index]) / t_r[index];

                for(unsigned int i = 0; i < 3; i++)
                {
                    if(i == index)
                        continue;

                    double valueC = p_r[i] + beta_2[0] * t_r[i];
                    double valueD = p_r[i] + beta_2[1] * t_r[i];

                    if(abs(valueC - C[i]) > tol || abs(valueD - D[i]) > tol)
                    {
                        sameLine = false;
                        break;
                    }
                }

                if(sameLine)
                {
                    std::sort(beta_1.begin(), beta_1.end());
                    std::sort(beta_2.begin(), beta_2.end());

                    return true;
                }

            }
        }

        return false;
    }
};

struct DFN{

    unsigned int NumberFractures = 0;
    std::vector<Fracture> Fractures = {};

    unsigned int NumberTraces = 0;
    std::vector<Trace> Traces = {};

};

bool FindIntersectionLine(const Eigen::Vector4d &plane1,
                          const Eigen::Vector4d &plane2,
                          Eigen::Vector3d &p_r,
                          Eigen::Vector3d &t_r,
                          const double &tol);

bool CompareTraces(const Trace &t1,
                   const Trace &t2);


bool ReadDFN(const std::string &fileName,
             DFN &dfn,
             const double &tol);

bool ImportFracture(const std::string &fileName,
                    DFN &dfn);

void CalculateTraces(DFN &dfn,
                     Fracture &f1,
                     Fracture &f2,
                     unsigned int &id,
                     const double &tol);

void WriteOutputFiles(const std::string &outputTracesFile,
                      const std::string &outputTipsFile,
                      const DFN &dfn);

}

#endif // DFN_HPP
