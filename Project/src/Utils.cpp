#include "fstream"
#include <iostream>
#include <sstream>
#include "Eigen/Eigen"
#include "iomanip"

#include "Utils.hpp"



bool ImportData(const std::string& inputFileName,
                unsigned int& nFracture,
                std::map<unsigned int, Eigen::MatrixXd>& Fractures)
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
            unsigned int FractureId;
            unsigned int NumVertices;
            std::istringstream converterId;
            std::istringstream convertN;
            std::getline(file, line, ';');
            converterId.str(line);
            converterId >> FractureId;
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

            Fractures[FractureId] = vertices; // ABBIAMO CREATO IL DIZIONARIO

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
    double majorDistance = 0;
    /*
    for (unsigned int i = 0; i < FirstFracture.cols(); i++)
    {
        double distance = firstCenter - FirstFracture.row(i);
        USARE LA NORMA 1... NO RADICE
        if ()
    }
    */
    return true;
}
