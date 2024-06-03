#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

#include "Utils.hpp"

int main()
{
    constexpr double tol = 10e-10;
    std::string inputFileName = "./DFN/FR3_data.txt";
    std::string outputFileName = "./DFN/TRACE_data.txt";
    std::string outputFileName2 = "./DFN/Output2.txt";

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
    unsigned int count = 0;
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
                    Data::Trace foundTrace;
                    if (FractureOperations::findTraces(Fractures[i], Fractures[j], t, foundTrace))
                    {                        
                        foundTrace.TraceId = (count)++;
                        Traces.push_back(foundTrace);
                        if (foundTrace.Tips[0] == false)
                            Fractures[i].passingTracesId.push_back(foundTrace.TraceId);
                        else
                            Fractures[i].notPassingTracesId.push_back(foundTrace.TraceId);
                        if (foundTrace.Tips[1] == false)
                            Fractures[j].passingTracesId.push_back(foundTrace.TraceId);
                        else
                            Fractures[j].notPassingTracesId.push_back(foundTrace.TraceId);

                    }

                }
            }
            else if (Fractures[i].normals == Fractures[j].normals && Fractures[i].d == Fractures[j].d ) // BOOK CASE
            {
                Data::Trace foundTrace;
                if (FractureOperations::bookCase(Fractures[i], Fractures[j], foundTrace))
                {
                    foundTrace.TraceId = (count)++;
                    Traces.push_back(foundTrace);
                    if (foundTrace.Tips[0] == false)
                        Fractures[i].passingTracesId.push_back(foundTrace.TraceId);
                    else
                        Fractures[i].notPassingTracesId.push_back(foundTrace.TraceId);
                    if (foundTrace.Tips[1] == false)
                        Fractures[j].passingTracesId.push_back(foundTrace.TraceId);
                    else
                        Fractures[j].notPassingTracesId.push_back(foundTrace.TraceId);
                }
                std::cout << " book case " << std::endl;
            }
        }
    }

    if (!Data::ExportData(outputFileName, Traces))
    {
        std::cerr<< "error: export failed"<< std::endl;
        return -2;
    }
    else
        std::cout << "Export successful " << std::endl;


    if (!Data::ExportSecondFile(outputFileName2, Fractures, Traces))
    {
        std::cerr<< "error: export failed"<< std::endl;
        return -3;
    }
    else
        std::cout << "Export successful " << std::endl;

    return 0;
}
