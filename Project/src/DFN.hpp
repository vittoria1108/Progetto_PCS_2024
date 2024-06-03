#ifndef DFN_HPP
#define DFN_HPP

#include <iostream>
#include <Eigen/Eigen>


using namespace std;
using namespace Eigen;


namespace FractureLibrary{

double CalculateSquareDistance(const Vector3d point1,
                               const Vector3d point2);

struct Trace{

    unsigned int Id;
    Vector2i FracturesIds = {};
    MatrixXd EndpointsCoordinates = {};
    double Length;

};

struct Fracture{

    unsigned int Id = 0;
    unsigned int NumberVertices = 0;
    MatrixXd VerticesCoordinates = {};
    Vector3d Barycentre = {};
    map<unsigned int, bool> Tips = {};

    vector<Trace> nTraces = {};
    vector<Trace> pTraces = {};

    Vector4d CalculatePlane()
    {
        Vector3d vector1 = VerticesCoordinates.col(1) - VerticesCoordinates.col(0);
        Vector3d vector2 = VerticesCoordinates.col(2) - VerticesCoordinates.col(0);

        Vector3d normal = vector1.cross(vector2);

        Vector3d point = VerticesCoordinates.col(0);
        double d = normal.dot(point);

        Vector4d plane = {normal[0], normal[1], normal[2], d};

        return plane;
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

    bool IntersectsLine(const Vector3d &p_r,
                        const Vector3d &t_r,
                        Vector2d &beta,
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

            Vector3d p_s = VerticesCoordinates.col(i);
            Vector3d t_s = VerticesCoordinates.col(next) - VerticesCoordinates.col(i);
            Vector3d prod = t_s.cross(t_r);

            if (prod.norm() > tol)
            {
                double alpha = (((p_r - p_s).cross(t_r)).dot(prod)) / (prod.dot(prod));

                if(alpha >= 0 && alpha < 1)
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
            else
            {
                bool sameLine = true;
                vector<unsigned int> indexBeta;
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
                    }
                }

                /* Se le coordinate dei punti sono uguali quando non riesco a calcolare le beta,
                    controllo che le beta calcolabili siano uguali*/
                if(sameLine)
                {
                    unsigned int numBeta = indexBeta.size() - 1;
                    oldBeta = (p_s[numBeta] - p_r[numBeta]) / t_r[numBeta];;

                    for(unsigned int j = 0; j < numBeta; j++)
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
                    beta[counter++] = betaTemp;

                    if(counter == 2)
                    {
                        sort(beta.begin(), beta.end());
                        return true;
                    }
                }
            }
        }

        return false;
    }
};

struct DFN{

    unsigned int NumberFractures = 0;
    vector<Fracture> Fractures = {};

    unsigned int NumberTraces = 0;
    vector<Trace> Traces = {};

};

bool FindIntersectionLine(const Vector4d &plane1,
                          const Vector4d &plane2,
                          Vector3d &p_r,
                          Vector3d &t_r,
                          const double &tol);

bool CompareTraces(const Trace &t1,
                   const Trace &t2);


bool ReadDFN(const string &fileName,
             DFN &dfn,
             const double &tol);

bool ImportFracture(const string &fileName,
                    DFN &dfn);

void CalculateTraces(DFN &dfn,
                     Fracture &f1,
                     Fracture &f2,
                     unsigned int &id,
                     const double &tol);

void WriteOutputFiles(const string &outputTracesFile,
                      const string &outputTipsFile,
                      const DFN &dfn);

}

#endif // DFN_HPP
