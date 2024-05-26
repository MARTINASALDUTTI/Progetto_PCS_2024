#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

#include "Utils.hpp"

int main()
{
    double tol = 10e-10;
    std::string inputFileName = "./DFN/FR3_data.txt";

    /*
     * creo la funzione import data a cui passo  nFracture e Fractures: un dizionario con chiave l'id e valore una struct
    contenente numero di vertici e coordinate
    */
    unsigned int nFracture = 0;
    //std::map<unsigned int, DFNlibrary::DFNdata> Fractures;
    std::vector<Data::Fract> Fractures;

    if (!Data::ImportData(inputFileName, nFracture, Fractures))
    {
        std::cerr<< "error: import failed"<< std::endl;
        return -1;
    }
    else
        std::cout << "Import successful" << std::endl;

    std::vector<Data::Trace> Traces;
    // Seleziono le fratture tra cui cercare le tracce e le trovo
    for (unsigned int i = 0; i < nFracture; i++)
    {
        for (unsigned int j = i+1; j < nFracture; j++)
        {
            //checking the distances between fractures
            if (FractureOperations::fracDistance(Fractures[i].vertices, Fractures[j].vertices))
                {
                    Eigen::Vector3d t = Fractures[i].normals.cross(Fractures[j].normals);
                    if (std::abs(t[0]) < tol && std::abs(t[1]) < tol && std::abs(t[2]) < tol) // checking if the two fractures are parallel
                    {
                        std::cout << "parallel" << std::endl;
                    }
                    else // else check book case
                    {
                        FractureOperations::findTraces(Fractures[i], Fractures[j], t);
                    }
                }
        }
    }
    return 0;
}