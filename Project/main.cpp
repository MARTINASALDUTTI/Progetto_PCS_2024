#include <iostream>

#include "Eigen/Eigen"

#include "Utils.hpp"

int main()
{
    constexpr double tol = 10e-10;
    std::string inputFileName = "./DFN/FR200_data.txt";
    std::string outputFileName = "./DFN/Trace_data.txt";
    std::string outputFileName2 = "./DFN/Trace_for_fracture.txt";

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

    std::list<Data::Trace> traces_list;
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

                if (t.squaredNorm() > tol)
                {
                    Data::Trace foundTrace;
                    if (FractureOperations::findTraces(Fractures[i], Fractures[j], t, foundTrace))
                    {
                        foundTrace.TraceId = (count)++;
                        traces_list.push_back(foundTrace);
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
                else if (t.squaredNorm() < tol) // BOOK CASE
                {
                    Data::Trace foundTrace;
                    if (FractureOperations::bookCase(Fractures[i], Fractures[j], foundTrace))
                    {
                        foundTrace.TraceId = (count)++;
                        traces_list.push_back(foundTrace);
                        if (foundTrace.Tips[0] == false)
                            Fractures[i].notPassingTracesId.push_back(foundTrace.TraceId);
                        else
                            Fractures[i].passingTracesId.push_back(foundTrace.TraceId);
                        if (foundTrace.Tips[1] == false)
                            Fractures[j].notPassingTracesId.push_back(foundTrace.TraceId);
                        else
                            Fractures[j].passingTracesId.push_back(foundTrace.TraceId);
                    }
                }
            }

        }
    }

    std::vector<Data::Trace> Traces;
    Traces.reserve(traces_list.size());
    for(Data::Trace& trace :traces_list)
    {
        Traces.push_back(trace);
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

    std::vector<PolygonalMeshLibrary::PolygonalMesh> Meshes;
    Meshes.reserve(Fractures.size());

    for (unsigned int i = 0; i < Fractures.size(); i++)
    {
        //per ogni frattura, salvo in un unico vettore tracce passanti ordinate per lunghezza decrescente +
        //tracce non ordinate per lunghezza decrescente

        PolygonalMeshLibrary::PolygonalMesh PolygonalMesh;
        //definisco ricorsivamente la funzione che fa i tagli
        if (Fractures[i].passingTracesId.size() + Fractures[i].notPassingTracesId.size()!= 0)
        {
            std::vector<Data::Trace> TracesCopy;
            TracesCopy = Traces;
            std::list<unsigned int>  AllTraces;
            AllTraces.insert(AllTraces.end(), Fractures[i].passingTracesId.begin(), Fractures[i].passingTracesId.end()); // Inserisci tutti gli elementi di vec1
            AllTraces.insert(AllTraces.end(), Fractures[i].notPassingTracesId.begin(), Fractures[i].notPassingTracesId.end());
            std::queue<Data::Fract> AllSubPolygons;
            AllSubPolygons.push(Fractures[i]);
            PolygonalMeshLibrary::MakeCuts(AllTraces,
                                           TracesCopy,
                                           PolygonalMesh,
                                           AllSubPolygons);
        }
        PolygonalMeshLibrary::CreateMesh(PolygonalMesh);

        PolygonalMesh.Num0DsCell = PolygonalMesh.coord0DsCellMap.size();
        PolygonalMesh.Num1DsCell = PolygonalMesh.Cell1DMap.size();
        PolygonalMesh.Num2DsCell = PolygonalMesh.Cell2DsVertices.size();

        Meshes.push_back(PolygonalMesh);
    }
    return 0;
}
