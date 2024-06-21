#include "fstream"
#include <iostream>
#include <sstream>
#include <limits>
#include "Eigen/Eigen"
#include <algorithm>
#include <iomanip>

#include "Utils.hpp"
#include "PolygonalMesh.hpp"

constexpr double tol=std::max(10*std::numeric_limits<double>::epsilon(), 10e-10);

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
        std::cerr<< "file open failed" << std::endl;
        return false;
    }
    Data::Fract fractData;

    //ignore line # Number of Fractures
    std::string line;
    std::getline(file, line);

    //read Number of Fractures
    std::getline(file, line);
    std::istringstream convert(line);
    convert >> nFracture;
    Fractures.reserve(nFracture);

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
            converterId.str(line);
            converterId >> fractData.FractId;
            std::getline(file,line);
            convertN.str(line);
            convertN >> NumVertices;

            //ignore # Vertices
            std::getline(file, line);

            //read the coordinates of the vertices
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

bool ExportData(const std::string& outputFileName,
                const std::vector<Data::Trace>& Traces)
{
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open())
    {
        std::cerr << "file open failed " << outputFileName << std::endl;
        return false;
    }

    outFile << "# Number of Traces\n";
    outFile << Traces.size() << "\n";

    // Scrivi i valori delle tracce
    for (unsigned int i = 0; i < Traces.size(); i++)
    {
        outFile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2\n";

        outFile << i << "; "
                << Traces[i].FractureIds[0] << "; "
                << Traces[i].FractureIds[1] << "; "
                << (Traces[i].ExtremesCoord[0])[0] << "; " << (Traces[i].ExtremesCoord[0])[1] << "; " << (Traces[i].ExtremesCoord[0])[2]<< "; "
                << (Traces[i].ExtremesCoord[1])[0] << "; " << (Traces[i].ExtremesCoord[1])[1] << "; " << (Traces[i].ExtremesCoord[1])[2] << "\n";
    }

    outFile.close();

    return true;
}

bool ExportSecondFile(const std::string& outputFileName,
                      std::vector<Data::Fract>& Fractures,
                      const std::vector<Data::Trace>& Traces)
{
    std::ofstream outFile(outputFileName);
    if (!outFile.is_open())
    {
        std::cerr << "file open failed " << outputFileName << std::endl;
        return false;
    }

    for (unsigned int i = 0; i < Fractures.size(); i++ )
    {
        if(Fractures[i].passingTracesId.size() + Fractures[i].notPassingTracesId.size() != 0)
        {
            outFile << "# FractureId; NumTraces\n";
            outFile << i << "; " << Fractures[i].passingTracesId.size() + Fractures[i].notPassingTracesId.size() << "\n";
            outFile << "# TraceId; Tips; Length\n";

            if (Fractures[i].passingTracesId.size() != 0)
            {
                SortLibrary::Mergesort(Fractures[i].passingTracesId, Traces);
                for(auto& elem :Fractures[i].passingTracesId )
                {
                    outFile << elem << "; " << 0 << "; " << Traces[elem].length << "\n";
                }
            }

            if( Fractures[i].notPassingTracesId.size() != 0)
            {
                SortLibrary::Mergesort(Fractures[i].notPassingTracesId, Traces);
                for(auto& elem :Fractures[i].notPassingTracesId )
                {
                    outFile << elem << "; " << 1 << "; " << Traces[elem].length << "\n";
                }
            }

        }
    }
    return true;
}
}

namespace FractureOperations{
bool fracDistance(const Eigen::MatrixXd &FirstFracture,
                  const Eigen::MatrixXd &SecondFracture)
{
    Eigen::Vector3d firstCenter = FirstFracture.rowwise().sum()/(FirstFracture.cols());
    Eigen::Vector3d secondCenter = SecondFracture.rowwise().sum()/(SecondFracture.cols());
    double majorDistance1 = 0 ;
    double majorDistance2 = 0;


    //vedere se fare una funzione
    for (unsigned int j = 0; j < FirstFracture.cols(); j++)
    {
        double distanceVC = (firstCenter - FirstFracture.col(j)).squaredNorm() ;
        if (distanceVC > majorDistance1)
            majorDistance1 = distanceVC;
    }

    for (unsigned int j = 0; j < SecondFracture.cols(); j++)
    {
        double distanceVC = (secondCenter - SecondFracture.col(j)).squaredNorm() ;
        if (distanceVC > majorDistance2)
            majorDistance2 = distanceVC;
    }

    double majorDistance = std::max(majorDistance1, majorDistance2);
    double centerDistance = (firstCenter - secondCenter).squaredNorm();

    if (centerDistance < 4*majorDistance)
        return true;

    return false;
}

void computePlane(Data::Fract& Fracture)
{
    Eigen::Vector3d vector1 = Fracture.vertices.col(0)-Fracture.vertices.col(1);
    Eigen::Vector3d vector2 =Fracture.vertices.col(0)-Fracture.vertices.col(2);

    Fracture.normals = (vector1.cross(vector2)).normalized();

    Fracture.d = Fracture.normals.dot(Fracture.vertices.col(0));
}


void findPosition(const Data::Fract& Fracture,
                  const Eigen::Vector3d& t,
                  const Eigen::Vector3d& P,
                  std::vector<Eigen::VectorXd>& CandidatePoints)
{
    double tol = 1e-10;
    bool previous = false;
    if (t.dot(Fracture.vertices.col(0) - P) >=  tol)
        previous = true;

    for (unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        bool current = false;
        if (t.dot(Fracture.vertices.col((i + 1) % Fracture.vertices.cols())-P)>=  tol)
            current = true;

        if(current != previous )
        {
            //if cosines have opposite sign, we have to solve the sistem
            Eigen::Vector3d V1 = Fracture.vertices.col(i);
            Eigen::Vector3d V2 = Fracture.vertices.col((i + 1) % Fracture.vertices.cols());

            Eigen::Vector3d intersection;
            if (FractureOperations::findExtreme(V1, V2,t, P, intersection))
            {
                CandidatePoints.push_back(intersection);
            }
        }
        if (((t.dot(Fracture.vertices.col((i + 1) % Fracture.vertices.cols())-P)>= -tol)||
             (t.dot(Fracture.vertices.col((i + 1) % Fracture.vertices.cols())-P)<= tol)) &&
            ((t.dot(Fracture.vertices.col(i)-P)>= -tol)||
             (t.dot(Fracture.vertices.col(i)-P))<= tol))
        {
            Eigen::Vector3d V1 = Fracture.vertices.col(i);
            Eigen::Vector3d V2 = Fracture.vertices.col((i + 1) % Fracture.vertices.cols());
            Eigen::Vector3d intersetion;
            if (FractureOperations::findExtreme(V1, V2,t, P, intersetion))
            {
                CandidatePoints.push_back(intersetion);
            }
        }
        previous = current;
    }
}

bool findExtreme(const Eigen::Vector3d& V1,
                 const Eigen::Vector3d& V2,
                 const Eigen::Vector3d& t,
                 const Eigen::Vector3d& P,
                 Eigen::Vector3d& intersection)
{
    double tol = 1e-10;
    Eigen::MatrixXd A(3,2);
    A.col(0) = t;
    A.col(1) = V1 - V2;

    Eigen::Vector3d b = V1 - P;

    Eigen::Vector2d paramVert = A.colPivHouseholderQr().solve(b);
    Eigen::Vector3d Candidate1= V1 + paramVert[1]*(V2-V1);
    Eigen::Vector3d Candidate2 = P + paramVert[0]* t;
    //intersection of the line t and V1-V2

    if ( ((Candidate1 - Candidate2).squaredNorm() < tol))
    {
        intersection = Candidate1;
        return true;
    }
    else
        return false;
}

bool isPointInPolygon(const Eigen::Vector3d& point,
                      const Eigen::MatrixXd& Polygon,
                      const Eigen::Vector3d& normal)
{
    bool PointInPolygon = true;
    for(unsigned int j = 0; j < Polygon.cols(); j++)
    {
        Eigen::Vector3d vector1 = Polygon.col((j+1) % Polygon.cols()) - Polygon.col(j);
        Eigen::Vector3d vector2 = point - Polygon.col(j);
        Eigen::Vector3d CrossProduct = vector1.cross(vector2);

        if(CrossProduct.dot(normal) < -tol )
            return !PointInPolygon;
    }

    return PointInPolygon;
}

bool isPointOnEdge(const Eigen::Vector3d& point, const Eigen::Vector3d& V1, const Eigen::Vector3d& V2)
{
    double tol = 10e-10;
    Eigen::Vector3d edge = V2 - V1;
    Eigen::Vector3d edgeToPoint = point - V1;
    double edgeLength = edge.norm();
    double projectedLength = edge.dot(edgeToPoint) / edgeLength;
    return (projectedLength >= 0 && projectedLength <= edgeLength && edgeToPoint.cross(edge).norm() < tol); 

}

bool isTracePassing(const Eigen::MatrixXd& fracture, const Eigen::Vector3d& traceStart, const Eigen::Vector3d& traceEnd)
{
    bool startOnEdge = false;
    bool endOnEdge = false;

    for (unsigned int i = 0; i < fracture.cols(); ++i) {
        Eigen::Vector3d V1 = fracture.col(i);
        Eigen::Vector3d V2 = fracture.col((i + 1) % fracture.cols());

        if (isPointOnEdge(traceStart, V1, V2)) {
            startOnEdge = true;
        }
        if (isPointOnEdge(traceEnd, V1, V2)) {
            endOnEdge = true;
        }
    }

    return startOnEdge && endOnEdge;
}

bool bookCase(const Data::Fract& FirstFracture,
              const Data::Fract& SecondFracture,
              Data::Trace& foundTrace)
{
    bool check = false;
    std::vector<Eigen::Vector3d> candidatePoints;
    for (unsigned int i = 0; i < FirstFracture.vertices.cols(); i++)
    {
        for(unsigned int j = 0; j < SecondFracture.vertices.cols(); j++)
        {
            if(isPointOnEdge(FirstFracture.vertices.col(i), SecondFracture.vertices.col(j), SecondFracture.vertices.col((j + 1) % SecondFracture.vertices.cols())))
            {
                candidatePoints.push_back(FirstFracture.vertices.col(i));
            }
        }
    }
    for (unsigned int i = 0; i < SecondFracture.vertices.cols(); i++)
    {
        for(unsigned int j = 0; j < FirstFracture.vertices.cols(); j++)
        {
            if(isPointOnEdge(SecondFracture.vertices.col(i), FirstFracture.vertices.col(j), FirstFracture.vertices.col((j + 1) % FirstFracture.vertices.cols())))
            {
                candidatePoints.push_back(SecondFracture.vertices.col(i));
            }
        }
    }
    if (candidatePoints.size() != 0)
    {
        std::vector<Eigen::VectorXd> extremePoints;
        extremePoints.push_back(candidatePoints[0]);

        for (unsigned int i = 1; i < candidatePoints.size(); i++)
        {
            if (candidatePoints[i] != candidatePoints[i-1])
                extremePoints.push_back(candidatePoints[i]);
        }

        if (extremePoints.size() == 2)
        {

        foundTrace.FractureIds[0] = FirstFracture.FractId;
        foundTrace.FractureIds[1] = SecondFracture.FractId;

        foundTrace.ExtremesCoord[0] = extremePoints[0];
        foundTrace.ExtremesCoord[1] = extremePoints[1];

        foundTrace.Tips[0] = !isTracePassing(FirstFracture.vertices, foundTrace.ExtremesCoord[0], foundTrace.ExtremesCoord[1]);
        foundTrace.Tips[1] = !isTracePassing(SecondFracture.vertices, foundTrace.ExtremesCoord[0], foundTrace.ExtremesCoord[1]);

        foundTrace.length = (foundTrace.ExtremesCoord[0] - foundTrace.ExtremesCoord[1]).norm();

        check = true;
        }
    }
    return check;
}

bool findTraces(const Data::Fract& FirstFracture,
                const Data::Fract& SecondFracture,
                const Eigen::Vector3d& t,
                Data::Trace& foundTrace)
{
    //find the intersection line -> equation r: p+tq

    Eigen::MatrixXd A(3, 3);
    A.row(0) = FirstFracture.normals;
    A.row(1) = SecondFracture.normals;
    A.row(2) = t;

    Eigen::Vector3d b;
    b.row(0) << FirstFracture.d;
    b.row(1) << SecondFracture.d;
    b.row(2) << 0;

    Eigen::Vector3d P = A.colPivHouseholderQr().solve(b);
    std::vector<Eigen::VectorXd> CandidatePoints;
    FractureOperations::findPosition(FirstFracture, t, P, CandidatePoints);
    FractureOperations::findPosition(SecondFracture, t, P, CandidatePoints);

    //std::cout <<" CandidatePoints.size() " << CandidatePoints.size() << std::endl;

    std::vector<Eigen::Vector3d> potentialPoints;

    for (unsigned int i = 0; i < CandidatePoints.size(); i++)
    {
        if (FractureOperations::isPointInPolygon(CandidatePoints[i], FirstFracture.vertices,FirstFracture.normals ) &&
            FractureOperations::isPointInPolygon(CandidatePoints[i], SecondFracture.vertices, SecondFracture.normals ))
        {
            potentialPoints.push_back(CandidatePoints[i]);
        }
    }

    if (potentialPoints.size() != 0)
    {
        std::vector<Eigen::VectorXd> extremePoints;
        extremePoints.push_back(potentialPoints[0]);

        for (unsigned int i = 1; i < potentialPoints.size(); i++)
        {
            if (potentialPoints[i] != potentialPoints[i-1])
                extremePoints.push_back(potentialPoints[i]);
        }

        foundTrace.FractureIds[0] = FirstFracture.FractId;
        foundTrace.FractureIds[1] = SecondFracture.FractId;

        foundTrace.ExtremesCoord[0] = extremePoints[0];
        foundTrace.ExtremesCoord[1] = extremePoints[1];

        foundTrace.Tips[0] = !isTracePassing(FirstFracture.vertices, foundTrace.ExtremesCoord[0], foundTrace.ExtremesCoord[1]);
        foundTrace.Tips[1] = !isTracePassing(SecondFracture.vertices, foundTrace.ExtremesCoord[0], foundTrace.ExtremesCoord[1]);

        foundTrace.length = (foundTrace.ExtremesCoord[0] - foundTrace.ExtremesCoord[1]).norm();

        return true;
    }
    else
        return false;
}
}

namespace SortLibrary
{
void merge(std::vector<unsigned int>& vecIdTraces, const std::vector<Data::Trace>& traces, size_t left, size_t center, size_t right)
{
    assert(right >= left);
    size_t i = left;
    size_t j = center+1;
    size_t k = 0;

    std::vector<unsigned int> tmp(right - left + 1);

    while (i <= center && j <= right) {
        if (traces[vecIdTraces[i]].length >= traces[vecIdTraces[j]].length)
        {
            assert(k < tmp.size());
            assert(i < vecIdTraces.size());
            tmp[k++] = vecIdTraces[i++];
        }
        else
        {
            assert(k < tmp.size());
            assert(j < vecIdTraces.size());
            tmp[k++] = vecIdTraces[j++];
        }
    }

    while (i <= center)
    {
        assert(k < tmp.size());
        assert(i < vecIdTraces.size());
        tmp[k++] = vecIdTraces[i++];
    }

    while (j <= right)
    {
        assert(k < tmp.size());
        assert(j < vecIdTraces.size());
        tmp[k++] = vecIdTraces[j++];
    }

    assert(k == (right - left + 1));

    for (size_t h = left; h <= right; h++)
    {
        vecIdTraces[h] = tmp[h-left];
    }
}

void mergesort(std::vector<unsigned int>& vecIdTraces, const std::vector<Data::Trace>& traces, size_t left, size_t right)
{
    assert(left <= vecIdTraces.size());
    assert(right <= vecIdTraces.size());

    if (left < right)
    {
        size_t center = (left + right)/2;
        mergesort(vecIdTraces, traces, left, center);
        mergesort(vecIdTraces, traces, center+1, right);
        merge(vecIdTraces, traces, left, center, right);
    }

}

void Mergesort(std::vector<unsigned int>& data, const std::vector<Data::Trace>& traces)
{
    if (data.size()>0)
    {
        mergesort(data, traces, 0, data.size()-1);
    }
}
}

namespace PolygonalMeshLibrary
{
bool SolveSystem(const Eigen::Vector3d& Direction,
                 const Eigen::Vector3d& vertice1,
                 const Eigen::Vector3d& vertice2,
                 const Eigen::Vector3d& point,
                 Eigen::Vector3d& Solution)
{
    bool Flag = false;

    Eigen::MatrixXd A(3,2);
    A.col(0) = Direction;
    A.col(1) = vertice1 - vertice2;

    Eigen::Vector3d b = vertice1 - point;

    Eigen::Vector2d paramVert = A.colPivHouseholderQr().solve(b);

    Eigen::Vector3d Candidate1= vertice1 + paramVert[1]*(vertice2-vertice1);
    Eigen::Vector3d Candidate2 = point + paramVert[0]* Direction;

    if ( ((Candidate1 - Candidate2).squaredNorm() < tol))
        {
            Solution = Candidate1;
            Flag = true;
        }

    return Flag;

}

bool MakeCuts(std::list<unsigned int>& AllTraces,
              std::vector<Data::Trace> &traces,
              PolygonalMeshLibrary::PolygonalMesh& PolygonMesh,
              std::queue<Data::Fract>& AllSubPolygons)
{
    Data::Fract CurrentPolygon;
    if (AllTraces.size() == 0)
    {
        while(!AllSubPolygons.empty())
        {
            PolygonalMeshLibrary::SavingSubpolygon(AllSubPolygons.front(), PolygonMesh);
            AllSubPolygons.pop();
        }
        return true;
    }
    else
    {
        if(AllSubPolygons.size() == 0 )
        {
            return false;
        }
        else
        {
            CurrentPolygon = AllSubPolygons.front();
        }

        //devo selezionare la traccia più lunga fra quelle del sottopoligono
        Data::Trace CurrentTrace;
        bool FineRicorsione = true;
        for (auto it = AllTraces.begin(); it != AllTraces.end(); it++)
        {
            if(FineRicorsione)
            {
                CurrentTrace = traces[*it];
                // vedere come gestire casi in cui la traccia viene divisa fra due sotto poligoni e chiedere se bisogna farlo
                if (FractureOperations::isPointInPolygon(CurrentTrace.ExtremesCoord[0], CurrentPolygon.vertices, CurrentPolygon.normals) ||
                    FractureOperations::isPointInPolygon(CurrentTrace.ExtremesCoord[1], CurrentPolygon.vertices, CurrentPolygon.normals))
                {
                        /*quando la traccia è interna al sottopoligono
                     * FineRicorsione diventa false ed esco dal ciclo
                     * altrimenti controllo se la traccia successiva è interns al poligonp
                     * Se non trovo questa traccia passo al poligonp suvccessivo
                     * FineRicorsione diventa true
                    */
                    bool TraceOnEdge = false;
                    FineRicorsione = checking(CurrentPolygon, CurrentTrace, TraceOnEdge);
                    if (TraceOnEdge)
                    {
                        auto punt = it;
                        it = AllTraces.erase(punt);
                    }
                }
                else if(PolygonalMeshLibrary::IsTraceInSubpolygon(CurrentPolygon, CurrentTrace, AllTraces, traces))
                {
                    bool TraceOnEdge = false;
                    FineRicorsione = checking(CurrentPolygon, CurrentTrace, TraceOnEdge);
                    if (TraceOnEdge)
                    {
                        auto punt = it;
                        it = AllTraces.erase(punt);
                    }
                }
            }
        }

        if(FineRicorsione == true)
        {
            // bisogna salvare il sotto poligono e uscire dalla funzione
            PolygonalMeshLibrary::SavingSubpolygon(CurrentPolygon,
                                                       PolygonMesh);

            // Rimuovi il primo sottopoligono analizzato dalla lista
            AllSubPolygons.pop();

            //richiamo makecut per ogni sotto pologono
            PolygonalMeshLibrary::MakeCuts(AllTraces,
                                           traces,
                                           PolygonMesh,
                                           AllSubPolygons);
        }
        else
        {
        //estraggo gli estremi e individuo la direzione sulla quale giace la traccia
        Eigen::Vector3d FirstExtreme = CurrentTrace.ExtremesCoord[0];
        Eigen::Vector3d SecondExtreme = CurrentTrace.ExtremesCoord[1];
        Eigen::Vector3d Direction = SecondExtreme - FirstExtreme;

        //per ogni vertice del poligono controllo la sua posizione reciproca rispetto alla traccia
        std::vector<Eigen::Vector3d> estremiTracce;
        estremiTracce.reserve(2);

        bool PreviousCheck = false;
        Eigen::Vector3d congiungente = CurrentPolygon.vertices.col(CurrentPolygon.vertices.cols()-1) - FirstExtreme;
        Eigen::Vector3d v= Direction.cross(congiungente);
        if (v.dot(CurrentPolygon.normals)>-tol)
            PreviousCheck = true;
        for (unsigned int i = 0; i < CurrentPolygon.vertices.cols(); i++)
        {
            /*per ogni lato ho tre opzioni
             * 1 il primo estremo è sul lato
             * 2 il secondo estremo è sul lato
             * 3 nessuno dei due estremi è sul lato, controllo se devo prolungare
             */
            bool CurrentCheck = false;

            Eigen::Vector3d congiungente = CurrentPolygon.vertices.col(i) - FirstExtreme;
            Eigen::Vector3d v= Direction.cross(congiungente);
            if (v.dot(CurrentPolygon.normals)>-tol)
                CurrentCheck = true;
            if(PreviousCheck != CurrentCheck)
            {
                if(FractureOperations::isPointOnEdge(FirstExtreme, CurrentPolygon.vertices.col(i), CurrentPolygon.vertices.col((i - 1) % CurrentPolygon.vertices.cols())))
                {
                    estremiTracce.push_back(FirstExtreme);
                }
                else if(FractureOperations::isPointOnEdge(SecondExtreme, CurrentPolygon.vertices.col(i), CurrentPolygon.vertices.col((i - 1) % CurrentPolygon.vertices.cols())))
                {
                    estremiTracce.push_back(SecondExtreme);
                }
                else if(PreviousCheck != CurrentCheck)
                {
                    unsigned int previusIndex;
                    if ( i == 0)
                    {
                        previusIndex =  CurrentPolygon.vertices.cols() - 1;
                    }
                    else
                    {
                        previusIndex = i-1;
                    }
                    //solve sistem
                    Eigen::Vector3d Solution;
                    if(PolygonalMeshLibrary::SolveSystem(Direction,
                                                          CurrentPolygon.vertices.col(i),
                                                          CurrentPolygon.vertices.col(previusIndex),
                                                          FirstExtreme,
                                                          Solution))
                        estremiTracce.push_back(Solution);
                }
            PreviousCheck = CurrentCheck;
            }
        }

        if(estremiTracce.size() < 2 )
        {
            if(estremiTracce.size() == 0)
            {
            PolygonalMeshLibrary::SavingSubpolygon(CurrentPolygon,
                                                   PolygonMesh);

            AllSubPolygons.pop();

            PolygonalMeshLibrary::MakeCuts(AllTraces,
                                           traces,
                                           PolygonMesh,
                                           AllSubPolygons);
            }

        }
        else
        {
            FirstExtreme = estremiTracce.front();
            Direction = estremiTracce.back() - FirstExtreme;
            std::vector<Eigen::Vector3d> FirstSide;
            std::vector<Eigen::Vector3d> SecondSide;

            unsigned int i = 0;
            congiungente = CurrentPolygon.vertices.col(i) - FirstExtreme;

            v= Direction.cross(congiungente);

            while (v.dot(CurrentPolygon.normals)> 0 )
            {
                if(v.dot(CurrentPolygon.normals)> 0)
                {
                    FirstSide.push_back(CurrentPolygon.vertices.col(i++));
                }
                else
                    SecondSide.push_back(CurrentPolygon.vertices.col(i++));

                congiungente = CurrentPolygon.vertices.col(i) - FirstExtreme;
                v= Direction.cross(congiungente);
            }
            //forse prima back e poi front
            FirstSide.push_back(estremiTracce.front());
            FirstSide.push_back(estremiTracce.back());
            SecondSide.push_back(estremiTracce.back());
            SecondSide.push_back(estremiTracce.front());

            for(unsigned int j = i; j < CurrentPolygon.vertices.cols(); j++)
            {
                Eigen::Vector3d congiungente = CurrentPolygon.vertices.col(j) - FirstExtreme;
                Eigen::Vector3d v= Direction.cross(congiungente);
                if (v.dot(CurrentPolygon.normals)> 0)
                {
                    FirstSide.push_back(CurrentPolygon.vertices.col(j));
                    //se la normale è parallela sta da un lato
                }
                else
                    SecondSide.push_back(CurrentPolygon.vertices.col(j));
            }

            //copio std::vector in Eigen::Matrix
            Eigen::MatrixXd SubPolygon(3,FirstSide.size());

            for(unsigned int k = 0; k < FirstSide.size(); k++)
            {
                SubPolygon.col(k) = FirstSide[k];
            }
            if(PolygonalMeshLibrary::CalculateArea(SubPolygon))
            {
                Data::Fract subpolygon;
                subpolygon.vertices = SubPolygon;
                //subpolygon.passingTracesId = CurrentPolygon.passingTracesId;
                subpolygon.normals = CurrentPolygon.normals;
                //aggiungo alla fine della lista i sottopoligoni da analizzare
                AllSubPolygons.push(subpolygon);
            }



            //copio std::vector in Eigen::Matrix
            Eigen::MatrixXd SecondSubPolygon(3,SecondSide.size());

            for(unsigned int k = 0; k < SecondSide.size(); k++)
            {
                SecondSubPolygon.col(k) = SecondSide[k];
            }

            if(PolygonalMeshLibrary::CalculateArea(SecondSubPolygon))
            {
                Data::Fract secondsubpolygon;
                secondsubpolygon.vertices = SecondSubPolygon;
                //secondsubpolygon.passingTracesId = CurrentPolygon.passingTracesId;
                secondsubpolygon.normals = CurrentPolygon.normals;

                //aggiungo alla fine della lista i sottopoligoni da analizzare
                AllSubPolygons.push(secondsubpolygon);
            }

            //elimino la traccia considerata:
            if (FractureOperations::isPointInPolygon(CurrentTrace.ExtremesCoord[0], CurrentPolygon.vertices, CurrentPolygon.normals) &&
                FractureOperations::isPointInPolygon(CurrentTrace.ExtremesCoord[1], CurrentPolygon.vertices, CurrentPolygon.normals))
            {
                AllTraces.remove(CurrentTrace.TraceId);
            }
            else if (FractureOperations::isPointInPolygon(CurrentTrace.ExtremesCoord[0], CurrentPolygon.vertices, CurrentPolygon.normals))
            {
                if((CurrentTrace.ExtremesCoord[1] - estremiTracce.back()).squaredNorm() < (CurrentTrace.ExtremesCoord[1] - estremiTracce.front()).squaredNorm() )
                    traces[CurrentTrace.TraceId].ExtremesCoord[0] = estremiTracce.back();
                else
                    traces[CurrentTrace.TraceId].ExtremesCoord[0] = estremiTracce.front();

                if ((CurrentTrace.ExtremesCoord[0] - CurrentTrace.ExtremesCoord[1]).squaredNorm() < tol)
                {
                    AllTraces.remove(CurrentTrace.TraceId);
                }
            }
            else if (FractureOperations::isPointInPolygon(CurrentTrace.ExtremesCoord[1], CurrentPolygon.vertices, CurrentPolygon.normals))
            {
                if((CurrentTrace.ExtremesCoord[0] - estremiTracce.back()).squaredNorm() < (CurrentTrace.ExtremesCoord[0] - estremiTracce.front()).squaredNorm() )
                {
                    traces[CurrentTrace.TraceId].ExtremesCoord[1] = estremiTracce.back();
                }
                else
                    traces[CurrentTrace.TraceId].ExtremesCoord[1] = estremiTracce.front();

                if ((CurrentTrace.ExtremesCoord[0] - CurrentTrace.ExtremesCoord[1]).squaredNorm() < tol)
                {
                    AllTraces.remove(CurrentTrace.TraceId);
                }
            }
            else
            {
                PolygonalMeshLibrary::UpdateTrace(CurrentTrace,
                                                  AllTraces,
                                                  traces,
                                                  estremiTracce);
            }

            traces[CurrentTrace.TraceId].length = (CurrentTrace.ExtremesCoord[0] - CurrentTrace.ExtremesCoord[0]).norm();

            // Rimuovi il sottopoligono analizzato dalla lista
            AllSubPolygons.pop();

            //richiamo makecut per ogni sotto pologono
            PolygonalMeshLibrary::MakeCuts(AllTraces,
                                           traces,
                                           PolygonMesh,
                                           AllSubPolygons);
        }
        }
    }

    return true;
}

void SavingSubpolygon(const Data::Fract& CurrentPolygon,
                      PolygonalMeshLibrary::PolygonalMesh& PolygonMesh)
{
    for (unsigned int i = 0; i < CurrentPolygon.vertices.cols(); i++)
    {
        PolygonMesh.Vertices_list.push_back(CurrentPolygon.vertices.col(i));
    }


    for (unsigned int i = 0; i < CurrentPolygon.vertices.cols(); i++)
    {
        std::array<Eigen::Vector3d, 2> edge = {CurrentPolygon.vertices.col(i), CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols())};
        PolygonMesh.edges_list.push_back(edge);
        std::array<Eigen::Vector3d, 2> edgeContr = {CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols()), CurrentPolygon.vertices.col(i)};
        PolygonMesh.edges_list.push_back(edgeContr);
    }

    PolygonMesh.subpolygons_list.push_back(CurrentPolygon.vertices);
}


void CreateMesh(PolygonalMeshLibrary::PolygonalMesh& PolygonMesh)
{
    auto itVertices_list = PolygonMesh.Vertices_list.begin();
    unsigned int IdCell0d = 0;
    while(itVertices_list != PolygonMesh.Vertices_list.end())
    {
        Eigen::Vector3d vertice = *itVertices_list;
        PolygonMesh.coord0DsCellMap.insert({IdCell0d++, vertice});
        PolygonMesh.Vertices_list.remove(*itVertices_list);
        itVertices_list = PolygonMesh.Vertices_list.begin();
    }

    auto itEdges_list = PolygonMesh.edges_list.begin();
    unsigned int IdCell1d = 0;
    while(itEdges_list != PolygonMesh.edges_list.end())
    {
        std::array<Eigen::Vector3d, 2> Edge = *(itEdges_list++);
        std::array<Eigen::Vector3d, 2> EdgeContr = *itEdges_list;
        Eigen::Vector3d FirstExtreme = Edge[0];
        Eigen::Vector3d SecondExtreme = Edge[1];

        auto itFirstExtreme = std::find_if(PolygonMesh.coord0DsCellMap.begin(), PolygonMesh.coord0DsCellMap.end(),
                               [&FirstExtreme](const std::pair<unsigned int, Eigen::Vector3d>& element) {
                                   return element.second == FirstExtreme;
                               });
        auto itSecondExtreme = std::find_if(PolygonMesh.coord0DsCellMap.begin(), PolygonMesh.coord0DsCellMap.end(),
                               [&SecondExtreme](const std::pair<unsigned int, Eigen::Vector3d>& element) {
                                   return element.second == SecondExtreme;
                               });
        unsigned int FirstId = itFirstExtreme->first;
        unsigned int SecondId = itSecondExtreme->first;
        std::array<unsigned int, 2> EdgeId = {FirstId, SecondId};
        PolygonMesh.Cell1DMap.insert({IdCell1d++, EdgeId});
        PolygonMesh.edges_list.remove(Edge);
        PolygonMesh.edges_list.remove(EdgeContr);

        itEdges_list = PolygonMesh.edges_list.begin();
    }

    unsigned int IdCell2d = 0;

    for(auto itVertices2D_list = PolygonMesh.subpolygons_list.begin(); itVertices2D_list != PolygonMesh.subpolygons_list.end(); itVertices2D_list++)
    {
        Eigen::MatrixXd matrix = *(itVertices2D_list);
        std::vector<unsigned int> Cell2DsVertices_vector;
        Cell2DsVertices_vector.reserve(matrix.cols());
        std::vector<unsigned int> Cell2DsEdges_vector;
        Cell2DsEdges_vector.reserve(matrix.cols() + 1);
        Eigen::Vector3d PreviousExtreme = matrix.col(matrix.cols() - 1);
        auto itPreviousExtreme = std::find_if(PolygonMesh.coord0DsCellMap.begin(), PolygonMesh.coord0DsCellMap.end(),
                                      [&PreviousExtreme](const std::pair<unsigned int, Eigen::Vector3d>& element) {
                                          return element.second == PreviousExtreme;
                                      });
        unsigned int PreviousId = itPreviousExtreme->first;

        for (unsigned int i = 0; i < matrix.cols(); i++ )
        {
            Eigen::Vector3d CurrentExtreme = matrix.col(i);
            auto itCurrentExtreme = std::find_if(PolygonMesh.coord0DsCellMap.begin(), PolygonMesh.coord0DsCellMap.end(),
                                               [&CurrentExtreme](const std::pair<unsigned int, Eigen::Vector3d>& element) {
                                                   return element.second == CurrentExtreme;
                                               });
            unsigned int CurrentId = itCurrentExtreme->first;
            Cell2DsVertices_vector.push_back(CurrentId);
            std::array<unsigned int, 2> Edge = { PreviousId, CurrentId };
            auto itEdge = std::find_if(PolygonMesh.Cell1DMap.begin(), PolygonMesh.Cell1DMap.end(),
                                       [&Edge](const std::pair<unsigned int, std::array<unsigned int, 2>>& element) {
                                     return element.second == Edge;
                                                 });
            if(itEdge == PolygonMesh.Cell1DMap.end())
            {
                std::array<unsigned int, 2> EdgeContr = { CurrentId, PreviousId };
                itEdge = std::find_if(PolygonMesh.Cell1DMap.begin(), PolygonMesh.Cell1DMap.end(),
                                     [&EdgeContr](const std::pair<unsigned int, std::array<unsigned int, 2>>& element) {
                                       return element.second == EdgeContr;
                                     });
            }
            unsigned int EdgeId = itEdge->first;
            Cell2DsEdges_vector.push_back(EdgeId);
            PreviousId = CurrentId;
        }

        PolygonMesh.Cell2DsVertices.insert({IdCell2d, Cell2DsVertices_vector});
        PolygonMesh.Cell2DsEdges.insert({IdCell2d, Cell2DsEdges_vector});
        IdCell2d++;
    }
}

bool checking(const Data::Fract& CurrentPolygon,
              Data::Trace& CurrentTrace,
              bool &TraceOnEdge)
{
    for(unsigned int i = 0; i < CurrentPolygon.vertices.cols(); i++)
    {
        if((CurrentPolygon.vertices.col(i) - CurrentTrace.ExtremesCoord[0]).squaredNorm() < tol ||
            (CurrentPolygon.vertices.col(i) - CurrentTrace.ExtremesCoord[1]).squaredNorm() < tol)
        {
            if(FractureOperations::isPointOnEdge(CurrentTrace.ExtremesCoord[0], CurrentPolygon.vertices.col(i), CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols())) &&
                FractureOperations::isPointOnEdge(CurrentTrace.ExtremesCoord[1], CurrentPolygon.vertices.col(i), CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols())) )
            {
                TraceOnEdge = true;
                return true;
            }
            else if (FractureOperations::isPointOnEdge(CurrentTrace.ExtremesCoord[0], CurrentPolygon.vertices.col(i), CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols())))
            {

                if((CurrentPolygon.vertices.col(i) - CurrentTrace.ExtremesCoord[1]).squaredNorm() < (CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols()) - CurrentTrace.ExtremesCoord[1]).squaredNorm())
                {
                    CurrentTrace.ExtremesCoord[0] = CurrentPolygon.vertices.col(i);
                    return true;
                }
                else if((CurrentPolygon.vertices.col(i) - CurrentTrace.ExtremesCoord[1]).squaredNorm() > (CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols()) - CurrentTrace.ExtremesCoord[1]).squaredNorm())
                {
                    CurrentTrace.ExtremesCoord[0] = CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols());
                    return true;
                }
            }
            else if (FractureOperations::isPointOnEdge(CurrentTrace.ExtremesCoord[1], CurrentPolygon.vertices.col(i), CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols())))
            {
                if((CurrentPolygon.vertices.col(i) - CurrentTrace.ExtremesCoord[0]).squaredNorm() < (CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols()) - CurrentTrace.ExtremesCoord[0]).squaredNorm())
                {
                    CurrentTrace.ExtremesCoord[1] = CurrentPolygon.vertices.col(i);
                    return true;
                }
                else if((CurrentPolygon.vertices.col(i) - CurrentTrace.ExtremesCoord[0]).squaredNorm() > (CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols()) - CurrentTrace.ExtremesCoord[0]).squaredNorm())
                {
                    CurrentTrace.ExtremesCoord[1] = CurrentPolygon.vertices.col((i + 1) % CurrentPolygon.vertices.cols());
                    return true;
                }
            }
        }
    }
    return false;
}

bool IsTraceInSubpolygon(const Data::Fract& CurrentPolygon,
                         Data::Trace& CurrentTrace,
                         std::list<unsigned int> &AllTraces,
                         std::vector<Data::Trace> &traces)
{
    bool PreviousCheck = false;
    Eigen::Vector3d Direction = CurrentTrace.ExtremesCoord[0] - CurrentTrace.ExtremesCoord[1];
    Eigen::Vector3d FirstExtreme = CurrentTrace.ExtremesCoord[1];
    Eigen::Vector3d vertice = CurrentPolygon.vertices.col(CurrentPolygon.vertices.cols() - 1) - FirstExtreme;
    Eigen::Vector3d CrossProduct = Direction.cross(vertice);
    if(CrossProduct.dot(CurrentPolygon.normals)> tol )
        PreviousCheck = true;

    for( unsigned int i = 0; i < CurrentPolygon.vertices.cols(); i++)
    {
        bool CurrentCheck = false;
        Eigen::Vector3d congiungente = CurrentPolygon.vertices.col(i) - CurrentTrace.ExtremesCoord[1];
        Eigen::Vector3d v= Direction.cross(congiungente);
        if (v.dot(CurrentPolygon.normals)>-tol)
            CurrentCheck = true;
        if(PreviousCheck != CurrentCheck)
        {
            //solve sistem
            Eigen::Vector3d Solution;
            if(PolygonalMeshLibrary::SolveSystem(Direction,
                                                 CurrentPolygon.vertices.col(i),
                                                 CurrentPolygon.vertices.col((i+1) % CurrentPolygon.vertices.cols()),
                                                 CurrentTrace.ExtremesCoord[1],
                                                 Solution))
            {
                if(FractureOperations::isPointOnEdge(Solution, CurrentPolygon.vertices.col(i), CurrentPolygon.vertices.col((i+1) % CurrentPolygon.vertices.cols())) &&
                   FractureOperations::isPointOnEdge(Solution, CurrentTrace.ExtremesCoord[0], CurrentTrace.ExtremesCoord[1]))
                {
                    return true;
                }
            }
        }
        PreviousCheck = CurrentCheck;
    }
    return false;
}

bool UpdateTrace(Data::Trace& CurrentTrace,
             std::list<unsigned int>& AllTraces,
             std::vector<Data::Trace>& traces,
             std::vector<Eigen::Vector3d> &estremiTracce)
{
    AllTraces.remove(CurrentTrace.TraceId);
    unsigned int max_id = traces.size();
    if((CurrentTrace.ExtremesCoord[0] - estremiTracce.front()).squaredNorm() < (CurrentTrace.ExtremesCoord[0] - estremiTracce.back()).squaredNorm())
    {
        Data::Trace FirstTraces;
        FirstTraces.TraceId = max_id;
        FirstTraces.ExtremesCoord[0] = CurrentTrace.ExtremesCoord[0] ;
        FirstTraces.ExtremesCoord[1] = estremiTracce.front();
        FirstTraces.length = (FirstTraces.ExtremesCoord[0] - FirstTraces.ExtremesCoord[1]).norm();
        if(FirstTraces.length > tol)
        {
            AllTraces.push_back(max_id);
            traces.push_back(FirstTraces);
        }
    }
    else
    {
        Data::Trace FirstTraces;
        FirstTraces.TraceId = max_id;
        FirstTraces.ExtremesCoord[0] = CurrentTrace.ExtremesCoord[0] ;
        FirstTraces.ExtremesCoord[1] = estremiTracce.back();
        FirstTraces.length = (FirstTraces.ExtremesCoord[0] - FirstTraces.ExtremesCoord[1]).norm();
        if(FirstTraces.length > tol)
        {
            AllTraces.push_back(max_id);
            traces.push_back(FirstTraces);
        }
    }

    if((CurrentTrace.ExtremesCoord[1] - estremiTracce.front()).squaredNorm() < (CurrentTrace.ExtremesCoord[1] - estremiTracce.back()).squaredNorm())
    {
        Data::Trace SecondTraces;
        SecondTraces.TraceId = max_id + 1;
        SecondTraces.ExtremesCoord[1] = CurrentTrace.ExtremesCoord[1] ;
        SecondTraces.ExtremesCoord[0] = estremiTracce.front();
        SecondTraces.length = (SecondTraces.ExtremesCoord[0] - SecondTraces.ExtremesCoord[1]).norm();
        if(SecondTraces.length > tol)
        {
            AllTraces.push_back(max_id + 1);
            traces.push_back(SecondTraces);
        }
    }
    else
    {
        Data::Trace SecondTraces;
        SecondTraces.TraceId = max_id + 1;
        SecondTraces.ExtremesCoord[1] = CurrentTrace.ExtremesCoord[1] ;
        SecondTraces.ExtremesCoord[0] = estremiTracce.back();
        SecondTraces.length = (SecondTraces.ExtremesCoord[0] - SecondTraces.ExtremesCoord[1]).norm();
        if(SecondTraces.length > tol)
        {
            AllTraces.push_back(max_id + 1);
            traces.push_back(SecondTraces);
        }
    }

    return true;
}

bool CalculateArea(const Eigen::MatrixXd& Cell)
{
    double SquaredArea = 0;
    Eigen::Vector3d refPointCoordinates = Cell.col(0);

    for(unsigned int i = 1; i < Cell.cols() - 1; i++)
    {
        Eigen::Vector3d firstEdge = Cell.col(i) - refPointCoordinates;
        Eigen::Vector3d secondEdge = Cell.col(i+1) - refPointCoordinates;

        double PartialArea = 0.25 * (firstEdge.cross(secondEdge)).squaredNorm();
        SquaredArea += PartialArea;
    }

    if (SquaredArea > std::pow(tol, 4))
    {
        return true;
    }

    return false;
}
}


