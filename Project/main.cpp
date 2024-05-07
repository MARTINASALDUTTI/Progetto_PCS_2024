#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

#include "Utils.hpp"


int main()
{
    std::string inputFileName = "FR3_data.txt";

    /*
     * creo la funzione import data a cui passo  nFracture e Fractures: un dizionario con chiave l'id e valore una struct
    contenente numero di vertici e coordinate
    */
    unsigned int nFracture = 0;
    //std::map<unsigned int, DFNlibrary::DFNdata> Fractures;
    std::map<unsigned int, Eigen::MatrixXd> Fractures;

    if (!ImportData(inputFileName, nFracture, Fractures))
    {
        std::cerr<< "error: import failed"<< std::endl;
        return -1;
    }
    /* CREARE FUNZIONE PER STAMPARE
    else
        std::cout << "Import successful: n= " << nFracture << ", Fractures= " << Fractures << std::endl;
    */

    return 0;
}