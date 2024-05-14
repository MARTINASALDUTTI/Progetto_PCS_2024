#include "fstream"
#include <iostream>
#include <sstream>
#include "Eigen/Eigen"
#include "iomanip"

#include "Utils.hpp"



bool ImportData(const std::string& inputFileName,
                unsigned int& nFracture,
                std::vector<Eigen::MatrixXd>& Fractures)
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
            Eigen::MatrixXd vertices(3, NumVertices);
            for (unsigned int i = 0; i < 3; i++)
            {
                std::getline(file,line);
                std::replace(line.begin(),line.end(),';',' ');
                std::istringstream convertCoord; //righe matrice vertici
                convertCoord.str(line);
                for (unsigned int j = 0; j < NumVertices; j++)
                    convertCoord >> vertices(i, j);

            }
            Fractures.push_back(vertices); // ABBIAMO CREATO IL VECTOR
        }
    }

    file.close();
    return true;
}


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

bool areParallel(Eigen::MatrixXd& FirstFracture,
                  Eigen::MatrixXd& SecondFracture)
{
    Eigen::Vector3d normal1;
    Eigen::Vector3d normal2;
    Eigen::Vector3d vector1 = FirstFracture.col(0)-FirstFracture.col(1);
    Eigen::Vector3d vector2 =FirstFracture.col(0)-FirstFracture.col(2);

    Eigen::Vector3d vector3 = SecondFracture.col(0)-SecondFracture.col(1);
    Eigen::Vector3d vector4 =SecondFracture.col(0)-SecondFracture.col(2);

    normal1 = (vector1.cross(vector2)).normalized();
    normal2 = (vector3.cross(vector4)).normalized();

    //vedi che tolleranza usare

    if(normal1.dot(normal2) == 1)
        return true;
    else
        return false;

    /*
    if(normal1.dot(normal2) < 1 + 1e-16 && normal1.dot(normal2)> 1- 1e16)
        return true;
    else
        return false;
*/
}
