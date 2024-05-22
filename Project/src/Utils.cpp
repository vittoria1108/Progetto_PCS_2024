#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "Utils.hpp"
#include "Eigen/Eigen"


using namespace std;
using namespace Eigen;


namespace FractureLibrary{

bool ImportFracture(const string &fileName, DFN &dfn)
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
    dfn.Fractures.resize(dfn.NumberFractures);

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        getline(file, header);
        //cout << header << endl;    // # FractureId; NumVertices

        getline(file, line, ';');
        istringstream idFractures(line);
        idFractures >> dfn.Fractures[i].Id ;

        getline(file, line);
        istringstream nVertices(line);
        nVertices >> dfn.Fractures[i].NumVertices;

        getline(file, header);
        //cout << header << endl;   // # Vertices
        dfn.Fractures[i].VerticesCoordinates.resize(3, dfn.Fractures[i].NumVertices);
        for(unsigned int j = 0; j < 3; j++)
        {
            getline(file, line);
            replace(line.begin(), line.end(), ';', ' ');
            //cout << line << endl;
            istringstream converter(line);

            for(unsigned int k = 0; k < dfn.Fractures[i].NumVertices; k++)
            {
                double coordinate;
                converter >> coordinate;
                dfn.Fractures[i].VerticesCoordinates(j, k) = coordinate;
            }

        }

        Vector3d barycentre = dfn.Fractures[i].VerticesCoordinates.rowwise().mean();
        dfn.Fractures[i].Barycentre = barycentre;
    }
    //cout << end

    //stampo id e coordinate per ogni frattura
    for(unsigned int a = 0; a < dfn.NumberFractures; a++)
    {
        cout << " id: " << dfn.Fractures[a].Id << endl;
        for(unsigned int b = 0; b < 3; b++)
        {
            cout << "[ ";
            for(unsigned int c = 0; c < dfn.Fractures[a].NumVertices; c++)
                cout <<fixed << setprecision(16)<< dfn.Fractures[a].VerticesCoordinates(b,c) << " ";
            cout << "]" << endl;
        }
    }

    file.close();
    return true;

}

Vector4d CalculatePlane(const Fracture& f)
{
    Vector3d vec1 = f.VerticesCoordinates.col(1) - f.VerticesCoordinates.col(0);
    //cout << "vec1  " << vec1[0] <<" "<< vec1[1] << " " << vec1[2] << endl;
    Vector3d vec2 = f.VerticesCoordinates.col(2) - f.VerticesCoordinates.col(0);
   // cout << "vec2  " << vec2[0] <<" "<< vec2[1] << " " << vec2[2] << endl;


    Vector3d normal = vec1.cross(vec2);
    double d = normal[0]*f.VerticesCoordinates(0,0) + normal[1]*f.VerticesCoordinates(1,0) + normal[2]*f.VerticesCoordinates(2,0);
    Vector4d coeff = {normal[0], normal[1], normal[2], d};
    //cout << "coeff plane  " << coeff[0] <<" "<< coeff[1] << " " << coeff[2] << " " << coeff[3] << endl;
    return coeff;
}

double CalculateDistance(const Vector3d point1, const Vector3d point2)
{
    double result = sqrt((point1-point2).transpose() * (point1 - point2));
    return result;
}

bool CalculateTraces(const Fracture &f1, const Fracture &f2)
{
    double tol = 10* numeric_limits<double>::epsilon();
    //cout << tol<< endl;
    double r1 = 0; //calcolo raggio della circonferenza che contiene il poligono
    for(unsigned int i = 0; i < f1.NumVertices; i++)
    {
        double newR = CalculateDistance(f1.Barycentre, f1.VerticesCoordinates.col(i));

        if(newR > r1)
            r1 = newR;
    }

    double r2 = 0;//calcolo raggio della circonferenza che contiene il poligono
    for(unsigned int i = 0; i < f2.NumVertices; i++)
    {
        double newR = CalculateDistance(f2.Barycentre, f2.VerticesCoordinates.col(i));

        if(newR > r2)
            r2 = newR;
    }
    double barDistances = CalculateDistance(f1.Barycentre, f2.Barycentre); //calcolo la distanza tra i due baricentri
    if(barDistances - (r1 + r2) > tol) //se la distanza tra i due baricentri è maggiore dei due raggi ritorno false
        return false;

    Vector4d plane1 = CalculatePlane(f1);
    Vector4d plane2 = CalculatePlane(f2);
    //cout << "plane1  " << plane1[0] <<" "<< plane1[1] << " " << plane1[2] <<" " << plane1[3] << endl;
    //cout << "plane2  " << plane2[0] <<" "<< plane2[1] << " " << plane2[2] <<" " << plane2[3] << endl;
    Vector3d normal1 = {plane1[0], plane1[1], plane1[2]};
    Vector3d normal2 = {plane2[0], plane2[1], plane2[2]};
    Vector3d t_r = normal1.cross(normal2);
    cout << "t_r  " << t_r[0] <<" "<< t_r[1] << " " << t_r[2] << endl;

    if(t_r.norm() < tol) //se t è vettore nullo
    {
        return false;
    }


    // Se non sono paralleli trovo la retta di intersezione
    Matrix3d normalMatrix;
    normalMatrix.row(0) = normal1;
    normalMatrix.row(1) = normal2;
    normalMatrix.row(2) = t_r;

    Vector3d constantTerms = {plane1[3], plane2[3], 0}; // uso piano con direzione t e passante per l'origine

    Vector3d p_r = normalMatrix.fullPivLu().solve(constantTerms); // sistema (normalMatrix)(xp yp zp) = constrantTerms
    cout << "p_r  " << p_r[0] <<" "<< p_r[1] << " " << p_r[2] << endl;

    // Controllo se la retta interseca le figure (interseca i piani, ma interseca le figure?)
    bool intersection = false; //potremmo sostituire con un contatore che conta i punti di intersezione trovati e si ferma a 2 (se abbiamo già i due punti non serve che giri sugli altri lati)

    Vector2d Beta_1 = {};
    unsigned int counter = 0;
    for(unsigned int i = 0; i < f1.NumVertices; i++)
    {
        unsigned int next; //per tenere conto anche dell'ultimo lato
        if (i == (f1.NumVertices -1))
            next = 0;
        else
            next = i+1;

    //segmento s: Ps + alfa*ts
    //retta    r: Pr + beta*tr
        Vector3d p_s = f1.VerticesCoordinates.col(i);
        //cout << "p_s  " << p_s[0] <<" "<< p_s[1] << " " << p_s[2] << endl;
        Vector3d t_s = f1.VerticesCoordinates.col(next) - f1.VerticesCoordinates.col(i);
        //cout << "t_s  " << t_s[0] <<" "<< t_s[1] << " " << t_s[2] << endl;
        Vector3d prod = t_s.cross(t_r);
        //cout << "prod  " << prod[0] <<" "<< prod[1] << " " << prod[2] << endl;

        if (prod.norm() > 0) //controllo che t_s e t_r non siano parallele/coincidenti, dovremo poi distinguere tra parallele e coincidenti e lavorare sulle coincidenti (se c'è appoggiato sopra un intero lato
        {
            //se non ho 0 cerco alfa

            double alfa = ((p_r - p_s).cross(t_r).dot(prod))/(prod.dot(prod));
            cout << "alfa1 m1 " << alfa << endl;

            if( alfa > 0 && alfa <= 1) //se c'è combinazione convessa vuol dire che il punto è interno al segmento e cerco alfa e beta
            {
                double beta = (((p_s - p_r).cross(t_s)).dot(-prod))/(prod.dot(prod));
                //cout << "beta1 m1 "<< beta << endl;
                Beta_1[counter] = beta;
                counter ++;

                if(counter == 2) //ho già trovato i due punti di intersezione
                    break;

                intersection = true;
            };
        }
    }

    if(!intersection)
        return false;



    Vector2d Beta_2 = {};
    counter = 0;

    for(unsigned int i = 0; i < f2.NumVertices; i++)
    {
        unsigned int next; //per tenere conto anche dell'ultimo lato
        if (i == (f2.NumVertices -1))
            next = 0;
        else
            next = i+1;

        Vector3d p_s = f2.VerticesCoordinates.col(i);
        //cout << "p_s  " << p_s[0] <<" "<< p_s[1] << " " << p_s[2] << endl;
        Vector3d t_s = f2.VerticesCoordinates.col(next) - f2.VerticesCoordinates.col(i);
        //cout << "t_s  " << t_s[0] <<" "<< t_s[1] << " " << t_s[2] << endl;
        Vector3d prod = t_s.cross(t_r);
        //cout << "prod  " << prod[0] <<" "<< prod[1] << " " << prod[2] << endl;

        if (prod.norm() > 0)
        {
            //se non ho 0 cerco alfa
            Vector3d diff = p_r - p_s;
            //cout << "pr - ps " << diff[0] <<" "<< diff[1] << " " << diff[2] << endl;
            Vector3d prod1 = (p_r - p_s).cross(t_r);
            //cout << "pr - ps x tr  " << prod1[0] <<" "<< prod1[1] << " " << prod1[2] << endl;
            //cout << "p_s  " << p_s[0] <<" "<< p_s[1] << " " << p_s[2] << endl;
            double alfa = (((p_r - p_s).cross(t_r)).dot(prod))/(prod.dot(prod));
            cout << "alfa2 m1 "<< alfa << endl;

            if( alfa > 0 && alfa < 1) //se c'è combinazione convessa vuol dire che il punto è interno al segmento e cerco alfa e beta
            {
                double beta = (((p_s - p_r).cross(t_s)).dot(-prod))/(prod.dot(prod));
                //cout << "beta2 m1 " << beta << endl;
                Beta_2[counter] = beta;
                counter ++;

                if(counter == 2) //ho già trovato i due punti di intersezione
                    break;

                intersection = true;
            }
        }

    }


    if(!intersection)
        return false;

    if(intersection)
    {
        cout << "Beta1 [ " ;
        cout << Beta_1[0] << ", " << Beta_1[1]<< " ]" << endl;
        cout << "Beta2 [ " ;
        cout << Beta_2[0] << ", " << Beta_2[1]<< " ]" << endl;
    }


    return true;
}

}



/*
    intersection = false;

    //opzione vecchia Ludo
    for(unsigned int i = 0; i < f2.NumVertices - 1; i++)
    {
        Vector3d dir = f2.VerticesCoordinates.col(i + 1) - f2.VerticesCoordinates.col(i);
        Vector3d b = point - f2.VerticesCoordinates.col(i);
        MatrixXd A;
        A.resize(3,2);
        A.col(0) = dir
        A.col(1) = - t_r

        Vector2d coeff = A.().solve(b);

        if(solution[0] >= 0 && solution[0] <= 1)
        {

        }

        }
    }

    if(!intersection)
        return false;



    //Se ho trovato intersezioni per entrambe le figure controllo che ci sia corrispondenza tra tali intersezioni
    //Controllo combinazioni lineari tra le due coppie di punti trovate sui lati delle figure
    //Vedi appunti, se un pt di int della frat 1 si trova tra i due pt della frat 2 allora sarà pt della traccia, e viceversa
    //Se entrambi i pt della frat 1 sono interni ai pt della frat 2 allora la traccia è passante per la frat 1, e viceversa
    //Dei quattro pt a nostra disposizione solo 2 sono i pt della traccia (se si ha una traccia)


    for(unsigned int i = 0; i<2; i++)
    {
        //prima forse controlliamo se alcuni pt delle due frat sono uguali, così nel caso evitiamo calcoli inutili

        //A, B vertici lato, P pt di inters controllo se AP = v*AB con v tra 0 e 1
        Vector3d AB = Pinter1[1] - Pinter1[0];
        Vector3d AP =  Pinter2[i] - Pinter1[0];
        // forma compatta if else     (condizione) ? valore_se_vero : valore_se_falso;
        double v_x = (AB[0] != 0) ? AP[0]/AB[0] : 0;
        double v_y = (AB[1] != 0) ? AP[1]/AB[1] : 0;
        double v_z = (AB[2] != 0) ? AP[2]/AB[2] : 0;
        double v = v_x!=0 ? v_x : (v_y!=0 ? v_y : v_z);
        if(v > 0 && v < 1)
        {
            ;
        }
    }
    */
