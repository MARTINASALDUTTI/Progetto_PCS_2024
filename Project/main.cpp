#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

#include "Utils.hpp"

int main()
{
    constexpr double tol = 10e-10;
    std::string inputFileName = "./DFN/FR10_data.txt";
    std::string outputFileName = "./DFN/TRACE_data.txt";
    std::string outputFileName2 = "./DFN/Output2.txt";

    /*
     * creo la funzione import data a cui passo  nFracture e Fractures: un dizionario con chiave l'id e valore una struct
    contenente numero di vertici e coordinate
    */
    unsigned int nFracture = 0;
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
                    //bisogna togliere la stampa
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
                        Fractures[i].notPassingTracesId.push_back(foundTrace.TraceId);
                    else
                        Fractures[i].passingTracesId.push_back(foundTrace.TraceId);
                    if (foundTrace.Tips[1] == false)
                        Fractures[j].notPassingTracesId.push_back(foundTrace.TraceId);
                    else
                        Fractures[j].passingTracesId.push_back(foundTrace.TraceId);
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
    //secondo me bisogna creare un std::vector di PolygonalMesh e riservare tanto spazio quanrto le fratture

    std::vector<PolygonalMeshLibrary::PolygonalMesh> Meshes;
    Meshes.reserve(Fractures.size());

    /*
    for (unsigned int i = 0; i < Fractures.size(); i++)
    {
        //per ogni frattura, salvo in un unico vettore tracce passanti ordinate per lunghezza decrescente +
        //tracce non ordinate per lunghezza decrescente

        PolygonalMeshLibrary::PolygonalMesh PolygonalMesh;
        //definisco ricorsivamente la funzione che fa i tagli
        std::cout << "frattura numero: " << i << std::endl;
        if (Fractures[i].passingTracesId.size() + Fractures[i].notPassingTracesId.size()!= 0)
        {
            std::list<unsigned int>  AllTraces;
            AllTraces.insert(AllTraces.end(), Fractures[i].passingTracesId.begin(), Fractures[i].passingTracesId.end()); // Inserisci tutti gli elementi di vec1
            AllTraces.insert(AllTraces.end(), Fractures[i].notPassingTracesId.begin(), Fractures[i].notPassingTracesId.end());
            std::queue<Data::Fract> AllSubPolygons;
            AllSubPolygons.push(Fractures[i]);
            PolygonalMeshLibrary::MakeCuts(AllTraces,
                                           Traces,
                                           PolygonalMesh,
                                           AllSubPolygons);
        }
        //in caso si può chiamare direttamente createmesh
        //oppure Cell0DId+1 (forse meglio);
        PolygonalMesh.Num0DsCell = PolygonalMesh.coord0DsCell.size();
        PolygonalMesh.Num1DsCell = PolygonalMesh.coord1DsCell.size();
        PolygonalMesh.Num2DsCell = PolygonalMesh.Cell2DsVertices.size();

        Meshes.push_back(PolygonalMesh);


        /*
        for (auto& elem : PolygonalMesh.coord0DsCell )
            std::cout << elem.transpose() << ", " << std::endl;
        std::cout << std::endl;

        unsigned int j = 0;
        for (auto& elem : PolygonalMesh.coord1DsCell_Id0DS )
            std::cout << j++ << " " << elem[0] << " " << elem[1] << std::endl;
        std::cout << std::endl;



        unsigned int j = 0;
        for (auto& polygon : PolygonalMesh.Cell2DsVertices )
        {
            std::cout << j++ << " ";
            for (auto& elem : polygon )
                std::cout << elem << ", ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
        */
    //}
    PolygonalMeshLibrary::PolygonalMesh PolygonalMesh;

    std::list<unsigned int>  AllTraces;

    AllTraces.insert(AllTraces.end(), Fractures[0].passingTracesId.begin(), Fractures[0].passingTracesId.end()); // Inserisci tutti gli elementi di vec1
    AllTraces.push_back(Fractures[0].notPassingTracesId[0]); // Inserisci tutti gli elementi di vec1
    AllTraces.push_back(Fractures[0].notPassingTracesId[1]); // Inserisci tutti gli elementi di vec1
    AllTraces.push_back(Fractures[0].notPassingTracesId[2]);
    AllTraces.push_back(Fractures[0].notPassingTracesId[3]);
    AllTraces.push_back(Fractures[0].notPassingTracesId[4]);
    AllTraces.push_back(Fractures[0].notPassingTracesId[5]);


    //AllTraces.insert(AllTraces.end(), Fractures[0].notPassingTracesId.begin(), Fractures[0].notPassingTracesId.end());
    std::queue<Data::Fract> AllSubPolygons;
    AllSubPolygons.push(Fractures[0]);
    PolygonalMeshLibrary::MakeCuts(AllTraces,
                                   Traces,
                                   PolygonalMesh,
                                   AllSubPolygons);



    for (auto& elem : PolygonalMesh.coord0DsCell )
        std::cout << elem.transpose() << ", " << std::endl;
    std::cout << std::endl;

    unsigned int j = 0;
    for (auto& elem : PolygonalMesh.coord1DsCell_Id0DS )
        std::cout << j++ << " " << elem[0] << " " << elem[1] << std::endl;
    std::cout << std::endl;


    //for (auto & elem : )

    /*
    unsigned int j = 0;
    for (auto& polygon : PolygonalMesh.Cell2DsVertices )
    {
        std::cout << j++ << " ";
        for (auto& elem : polygon )
            std::cout << elem << ", ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
    */

    return 0;
}
