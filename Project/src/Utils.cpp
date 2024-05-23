#include "fstream"
#include <iostream>
#include <sstream>
#include "Eigen/Eigen"
#include "iomanip"

#include "Utils.hpp"


namespace Data{
bool ImportData(const std::string& inputFileName,
                unsigned int& nFracture,
                std::vector<Data::Fract>& Fractures)
{
    //Open file
    std::ifstream file;
    file.open(inputFileName);

    if (file.fail())
    {
        std::cerr<< "file open failed"<< std::endl;
        return false;
    }

    //ignore line # Number of Fractures
    std::string line;
    std::getline(file, line);

    //read Number of Fractures
    std::getline(file, line);
    std::istringstream convert(line);
    convert >> nFracture;

    while (!file.eof())
    {
        //ignore # FractureId; NumVertices
        std::getline(file, line);

        if (line != "")
        {

            //read FractureId NumVertices
            unsigned int NumVertices;
            std::istringstream converterId;
            std::istringstream convertN;
            std::getline(file, line, ';');
            std::getline(file,line);
            convertN.str(line);
            convertN >> NumVertices;

            //ignore # Vertices
            std::getline(file, line);

            //read the coordinates of the vertices
            Data::Fract fractData;
            Eigen::MatrixXd v(3, NumVertices); //v is the matrix of vertices' coordinates
            for (unsigned int i = 0; i < 3; i++)
            {
                std::getline(file,line);
                std::replace(line.begin(),line.end(),';',' ');
                std::istringstream convertCoord; //righe matrice vertici
                convertCoord.str(line);
                for (unsigned int j = 0; j < NumVertices; j++)
                    convertCoord >> v(i, j);

            }
            fractData.vertices = v;
            FractureOperations::computePlane(fractData);
            Fractures.push_back(fractData); // ABBIAMO CREATO IL VECTOR
        }
    }

    file.close();
    return true;
}
}

namespace FractureOperations{
bool fracDistance(Eigen::MatrixXd& FirstFracture,
                  Eigen::MatrixXd& SecondFracture)
{
    Eigen::Vector3d firstCenter = FirstFracture.rowwise().sum()/(FirstFracture.cols());
    Eigen::Vector3d secondCenter = SecondFracture.rowwise().sum()/(SecondFracture.cols());
    double majorDistance1 = 0;
    double majorDistance2 = 0;

    //vedere se fare una funzione
    for (unsigned int j = 0; j < FirstFracture.cols(); j++)
    {
        double  distanceVC = 0;
        for (unsigned int i = 0; i < 3; i++)
        {
            distanceVC += ((FirstFracture(i,j)- firstCenter(i)) * (FirstFracture(i,j)- firstCenter(i)));
        }
        if (distanceVC > majorDistance1)
            majorDistance1 = distanceVC;
    }

    for (unsigned int j = 0; j < SecondFracture.cols(); j++)
    {
        double  distanceVC = 0;
        for (unsigned int i = 0; i < 3; i++)
        {
            distanceVC += ((SecondFracture(i,j)- secondCenter(i)) * (SecondFracture(i,j)- secondCenter(i)));
        }
        if (distanceVC > majorDistance2)
            majorDistance2 = distanceVC;
    }

    double centerDistance = 0;
    for (unsigned int i = 0; i < 3; i++)
        centerDistance += (firstCenter(i) - secondCenter(i)) * (firstCenter(i) - secondCenter(i));
    if (centerDistance > majorDistance1 + majorDistance2)
        return false;
    else
        return true;
}

void computePlane(Data::Fract& Fracture)
{
    Eigen::Vector3d vector1 = Fracture.vertices.col(0)-Fracture.vertices.col(1);
    Eigen::Vector3d vector2 =Fracture.vertices.col(0)-Fracture.vertices.col(2);

    Fracture.normals = (vector1.cross(vector2)).normalized();

    Fracture.d = Fracture.normals.dot(Fracture.vertices.col(0));

}

void findTraces(const Data::Fract& FirstFracture,
                const Data::Fract& SecondFracture,
                const Eigen::Vector3d& t,
                std::vector<Data::Trace>& Traces)
{
    //chiedi  se escludere caso in cui l'intersezione è un solo punto
    //find the intersection line -> equation r: p+tq

    Eigen::MatrixXd A(3, 3);
    A.row(0) = FirstFracture.normals;
    A.row(1) = SecondFracture.normals;
    A.row(2) = t;

    Eigen::Vector3d b;
    b.row(0) << FirstFracture.d;
    b.row(1) << SecondFracture.d;
    b.row(2) << 0;

    Eigen::Vector3d P = A.colPivHouseholderQr().solve(b);

    std::cout << "Normale 1: " << FirstFracture.normals.transpose() << std::endl;
    std::cout << "Normale 2: " << SecondFracture.normals.transpose() << std::endl;
    std::cout << "Vettore t: " << t.transpose() << std::endl;
    std::cout << "Punto P: " << P.transpose() << std::endl;
    std::cout << std::endl;


    Eigen::Vector3d ExtremeTrace;

    std::cout << "FirstFracture" << std::endl;
    FractureOperations::findPosition(FirstFracture, t, P, ExtremeTrace);
    std::cout << "SecondFracture" << std::endl;

    FractureOperations::findPosition(SecondFracture, t, P, ExtremeTrace);

    }

void findPosition(const Data::Fract& FirstFracture,
                  const Eigen::Vector3d& t,
                  const Eigen::Vector3d& P,
                  Eigen::Vector3d& ExtremeTrace)
{

    double tol = 1e-10;
    bool previous = false;
    if (t.dot(FirstFracture.vertices.col(0)-P)>=- tol)
        previous = true;

    for (unsigned int i = 0; i < FirstFracture.vertices.cols(); i++)
    {
        bool current = false;
        if (t.dot(FirstFracture.vertices.col(i)-P)>= - tol)
            current = true;

        if(current != previous )
        {
            //if cosines have opposite sign, we have to solve the sistem
            Eigen::Vector3d V1 = FirstFracture.vertices.col(i);
            Eigen::Vector3d V2 = FirstFracture.vertices.col((i - 1) % FirstFracture.vertices.cols()); //-1%N = N-1 where N is a natural number
            FractureOperations::findExtreme(V1, V2,t, P, ExtremeTrace);
        }

        if (((t.dot(FirstFracture.vertices.col(i)-P)>= -tol)|| (t.dot(FirstFracture.vertices.col(i)-P)<= tol)) &&
         ((t.dot(FirstFracture.vertices.col((i - 1) % FirstFracture.vertices.cols())-P)>= -tol)||
          (t.dot(FirstFracture.vertices.col((i - 1) % FirstFracture.vertices.cols())-P))<= tol))
        {
            Eigen::Vector3d V1 = FirstFracture.vertices.col(i);
            Eigen::Vector3d V2 = FirstFracture.vertices.col((i - 1) % FirstFracture.vertices.cols());
            FractureOperations::findExtreme(V1, V2,t, P, ExtremeTrace);
        }

        previous = current;
    }
}

void findExtreme(const Eigen::Vector3d& V1,
                 const Eigen::Vector3d& V2,
                 const Eigen::Vector3d& t,
                 const Eigen::Vector3d& P,
                 Eigen::Vector3d& ExtremeTrace)
{
    double tol = 1e-10;
    Eigen::MatrixXd A(3,2);
    A.col(0) = t;
    A.col(1) = V1 - V2;

    Eigen::Vector3d b;
    b.row(0) << V1[0]-P[0];
    b.row(1) << V1[1]-P[1];
    b.row(2) << V1[2]-P[2];

    Eigen::Vector2d paramVert = A.colPivHouseholderQr().solve(b);
    Eigen::Vector3d Candidate1= V1 + paramVert[1]*(V2-V1);
    Eigen::Vector3d Candidate2 = P + paramVert[0]* t;

    if ( ((Candidate1 - Candidate2)[0] < tol && (Candidate1 - Candidate2)[0] > -tol) &&
         ((Candidate1 - Candidate2)[1] < tol && (Candidate1 - Candidate2)[1] > -tol) &&
         ((Candidate1 - Candidate2)[2] < tol && (Candidate1 - Candidate2)[2] > -tol))
    {
        std::cout << "paramVert " << paramVert.transpose() << std::endl;
        ExtremeTrace = Candidate1;
        std::cout << "ExtremeTrace " << ExtremeTrace.transpose() << std::endl;
        std::cout << std::endl;
    }
}
}


