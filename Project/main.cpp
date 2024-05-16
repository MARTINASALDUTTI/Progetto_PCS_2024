#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

#include "Utils.hpp"



int main()
{
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
            if (FractureOperations::fracDistance(Fractures[i].vertices, Fractures[j].vertices))
            {
                //checking if the fractures are parallel
                //vedi che tolleranza usare

                if(Fractures[i].normals.dot(Fractures[j].normals) != 1)
                {

                    FractureOperations::findTraces(Fractures[i].vertices, Fractures[j].vertices, Traces);
                    std::cout << "è andato zi" << std::endl;

                }
            }
        }
    }
    return 0;
}
