#include "fstream"
#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

#include "Utils.hpp"



bool ImportData(const std::string inputFileName,
                unsigned int nFracture,
                std::map<unsigned int, Eigen::MatrixXd> Fractures)
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
    std::istringstream convert(line);
    convert >> nFracture;

    /*
    //salvo tutte le righe per sostituire ; con " ". per le prime tre righe no perchè non ci sono  ; oppure ignoriamo
    std::list<std::string> listLines;
    std::string line;
    while (std::getline(file, line))
    {
        std::replace(line.begin(),line.end(),';',' ');
        listLines.push_back(line);
    }
    file.close();

    for (const std::string& line : listLines)
    {

    }
    */
    while (!file.eof())
    {
        //ignore # FractureId; NumVertices
        std::getline(file, line);

        //read FractureId NumVertices
        unsigned int FractureId;
        unsigned int NumVertices;
        std::getline(file, line, ';');
        std::istringstream converterId(line);
        converterId >> FractureId;
        std::istringstream convertN(line.substr(2));
        convertN >> NumVertices;

        std::cout << FractureId << " " << NumVertices<< std::endl;
        /*
        auto ret = mesh.Cell0DMarkers.insert({marker, {id}});
        if(!ret.second)
            (ret.first)->second.push_back(id);
        */
    }

    return true;
}