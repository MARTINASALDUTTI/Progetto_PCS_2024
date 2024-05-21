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

    Fracture.d = -Fracture.normals.dot(Fracture.vertices.col(0));

}

void findTraces(const Data::Fract& FirstFracture,
                const Data::Fract& SecondFracture,
                std::vector<Data::Trace>& Traces)
{
    //chiedi  se escludere caso in cui l'intersezione è un solo punto
    //find the intersection line -> equation r: p+tq

    Eigen::Vector3d t = FirstFracture.normals.cross(SecondFracture.normals);

    Eigen::MatrixXd A(3, 3);
    A.row(0) = FirstFracture.normals;
    A.row(1) = SecondFracture.normals;
    A.row(2) = t;

    Eigen::Vector3d b;
    b.row(0) << FirstFracture.d;
    b.row(1) << SecondFracture.d;
    b.row(2) << 0;

    Eigen::Vector3d P = A.colPivHouseholderQr().solve(b);

    //checking if the vertice and the intersection line make an angle with cos > 0 (true) else false
    bool previous = false;
    if (t.dot(FirstFracture.vertices.col(0)-P)>0)
        previous = true;

    //itero fino a FirstFracture.vertices.cols-1 perchè devo copnfrontare l'ultimo vertice con il primo e non con i+1
    for (unsigned int i = 0; i < FirstFracture.vertices.cols(); i++)
    {
        bool current = false;
        if (t.dot(FirstFracture.vertices.col(i)-P)>0 )
            current = true;

        if(current != previous)
        {
            //if cosines have opposite sign, we have to solve the sistem
            Eigen::Vector3d V1 = FirstFracture.vertices.col(i);
            Eigen::Vector3d V2 = FirstFracture.vertices.col((i + 1) % FirstFracture.vertices.cols());

            // eq retta:
            // eq piano:
            double t=-(SecondFracture.normals.dot(V1)+ SecondFracture.d) /(SecondFracture.normals.dot(V2 - V1));

            Eigen::Vector3d intersection = V1+t*(V2-V1);

            std::cout << intersection << std::endl;

            /*
            Eigen::Vector3d edge = V2 - V1;
            Eigen::Matrix2d A2d;
            A2d << t.head<2>(), -edge.head<2>();
            Eigen::Vector2d b2d = V1.head<2>() - P.head<2>();

            // Risolvi il sistema lineare per trovare i parametri t e s
            if (A2d.determinant() != 0)
            {
                Eigen::Vector2d ts = A2d.colPivHouseholderQr().solve(b2d);
                std::cout <<  ts  << std::endl;
                std::cout << std::endl;

                double t_param = ts(0);
                double s_param = ts(1);

                // Controlla se il punto di intersezione è all'interno del segmento
                if (s_param >= 0 && s_param <= 1)
                {
                    std::cout << " pippo " << std::endl;
                    Eigen::Vector3d intersection = P + t_param * t;
                        //intersectionPoints.push_back(intersection);
                        std::cout << intersection << std::endl;

                }

            }
*/
        }
        previous = current;
    }


}
}
