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
    Data::Fract fractData;

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
            converterId.str(line);
            converterId >> fractData.FractId;
            std::getline(file,line);
            convertN.str(line);
            convertN >> NumVertices;

            //ignore # Vertices
            std::getline(file, line);

            //read the coordinates of the vertices
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
                const Eigen::Vector3d& t)
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
    std::vector<Eigen::VectorXd> CandidatePoints;
    //Eigen::Vector3d ExtremeTrace;
    FractureOperations::findPosition(FirstFracture, t, P, CandidatePoints);
    FractureOperations::findPosition(SecondFracture, t, P, CandidatePoints);

    std::vector<Eigen::VectorXd> extremePoints;
    for (unsigned int i = 0; i < CandidatePoints.size(); i++)
    {
        if (isPointInPolygon(CandidatePoints[i], FirstFracture.vertices) && isPointInPolygon(CandidatePoints[i], SecondFracture.vertices))
            extremePoints.push_back(CandidatePoints[i]);
    }

    for (auto& elem : extremePoints)
        std::cout << elem.transpose() << std::endl;
    std::cout << std::endl;
    }

void findPosition(const Data::Fract& Fracture,
                  const Eigen::Vector3d& t,
                  const Eigen::Vector3d& P,
                  std::vector<Eigen::VectorXd>& CandidatePoints)
{
    double tol = 1e-10;
    bool previous = false;
    if (t.dot(Fracture.vertices.col(0)-P)>=- tol)
        previous = true;

    for (unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        bool current = false;
        if (t.dot(Fracture.vertices.col(i)-P)>= - tol)
            current = true;

        if(current != previous )
        {
            //if cosines have opposite sign, we have to solve the sistem
            Eigen::Vector3d V1 = Fracture.vertices.col(i);
            Eigen::Vector3d V2 = Fracture.vertices.col((i - 1) % Fracture.vertices.cols()); //-1%N = N-1 where N is a natural number
            Eigen::Vector3d intersetion;
            if (FractureOperations::findExtreme(V1, V2,t, P, intersetion))
            {
                CandidatePoints.push_back(intersetion);
            }
        }

        if (((t.dot(Fracture.vertices.col(i)-P)>= -tol)|| (t.dot(Fracture.vertices.col(i)-P)<= tol)) &&
         ((t.dot(Fracture.vertices.col((i - 1) % Fracture.vertices.cols())-P)>= -tol)||
          (t.dot(Fracture.vertices.col((i - 1) % Fracture.vertices.cols())-P))<= tol))
        {
            Eigen::Vector3d V1 = Fracture.vertices.col(i);
            Eigen::Vector3d V2 = Fracture.vertices.col((i - 1) % Fracture.vertices.cols());
            Eigen::Vector3d intersetion;
            if (FractureOperations::findExtreme(V1, V2,t, P, intersetion))
            {
                CandidatePoints.push_back(intersetion);
            }
        }
        previous = current;
    }
}

bool findExtreme(const Eigen::Vector3d& V1,
                 const Eigen::Vector3d& V2,
                 const Eigen::Vector3d& t,
                 const Eigen::Vector3d& P,
                 Eigen::Vector3d& intersection)
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
    //intersection of the line t and V1-V2

    if ( ((Candidate1 - Candidate2)[0] < tol && (Candidate1 - Candidate2)[0] > -tol) &&
         ((Candidate1 - Candidate2)[1] < tol && (Candidate1 - Candidate2)[1] > -tol) &&
         ((Candidate1 - Candidate2)[2] < tol && (Candidate1 - Candidate2)[2] > -tol))
    {
        intersection = Candidate1;
        return true;
    }
    else
        return false;
}


bool isPointInPolygon(const Eigen::Vector3d& point,
                      const Eigen::MatrixXd& Fracture)
{
    double tol = 1e-10;

    Eigen::Vector3d AB = Fracture.col(1) - Fracture.col(0);
    Eigen::Vector3d AD = Fracture.col(3) - Fracture.col(0);
    Eigen::Vector3d AP = point - Fracture.col(0);

    double dotAB_AB = AB.dot(AB);
    double dotAD_AD = AD.dot(AD);
    double dotAP_AB = AP.dot(AB);
    double dotAP_AD = AP.dot(AD);

    double lambda1 = dotAP_AB / dotAB_AB;
    double lambda2 = dotAP_AD / dotAD_AD;

    return (lambda1 >= -tol && lambda1 <= 1.0+tol && lambda2 >= -tol && lambda2 <= 1+tol);
    /*
    bool control = true;
    double tol = 1e-10;
    Eigen::VectorXd lambda;

    Eigen::MatrixXd A(4,4);

    A.row(0) = Fracture.row(0);
    A.row(1) = Fracture.row(1);
    A.row(2) = Fracture.row(2);
    A.row(3) << 1, 1, 1, 1;

    std::cout << A << std::endl;

    Eigen::Vector4d b;
    b << point(0), point(1), point(2), 1;
    std::cout << b << std::endl;

    b.row(0) = point.row(0);
    b.row(1) = point.row(1);
    b.row(2) = point.row(2);
    b.row(3) << 1;

    lambda = A.colPivHouseholderQr().solve(b);

    std::cout << lambda.transpose() << std::endl;

    for(unsigned int i = 0; i < lambda.rows(); i++)
    {
        if (lambda(i) < -tol)
            control = false;
    }
    std::cout << "control" << control << std::endl;
    return control;
*/

}
}



