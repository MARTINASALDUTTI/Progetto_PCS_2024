#include "fstream"
#include <iostream>
#include <sstream>
#include "Eigen/Eigen"
#include "iomanip"
#include <algorithm> // Per std::find
#include <unordered_map>


#include "Utils.hpp"
#include "PolygonalMesh.hpp"

constexpr double tol = 10e-10;

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
    if (t.dot(Fracture.vertices.col(0)-P)>=- tol)
        previous = true;

    for (unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        bool current = false;
        if (t.dot(Fracture.vertices.col(i)-P)>= - tol)
            current = true;

        if(current != previous )
        {
            //if cosines have opposite sign, we have to solve the sistem
            Eigen::Vector3d V1 = Fracture.vertices.col(i);
            Eigen::Vector3d V2 = Fracture.vertices.col((i - 1) % Fracture.vertices.cols()); //-1%N = N-1 where N is a natural number
            Eigen::Vector3d intersetion;
            if (FractureOperations::findExtreme(V1, V2,t, P, intersetion))
            {
                CandidatePoints.push_back(intersetion);
            }
        }

        if (((t.dot(Fracture.vertices.col(i)-P)>= -tol)|| (t.dot(Fracture.vertices.col(i)-P)<= tol)) &&
            ((t.dot(Fracture.vertices.col((i - 1) % Fracture.vertices.cols())-P)>= -tol)||
             (t.dot(Fracture.vertices.col((i - 1) % Fracture.vertices.cols())-P))<= tol))
        {
            Eigen::Vector3d V1 = Fracture.vertices.col(i);
            Eigen::Vector3d V2 = Fracture.vertices.col((i - 1) % Fracture.vertices.cols());
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

    Eigen::Vector3d b;
    b.row(0) << V1[0]-P[0];
    b.row(1) << V1[1]-P[1];
    b.row(2) << V1[2]-P[2];

    Eigen::Vector2d paramVert = A.colPivHouseholderQr().solve(b);
    Eigen::Vector3d Candidate1= V1 + paramVert[1]*(V2-V1);
    Eigen::Vector3d Candidate2 = P + paramVert[0]* t;
    //intersection of the line t and V1-V2

    if ( ((Candidate1 - Candidate2)[0] < tol && (Candidate1 - Candidate2)[0] > -tol) &&
        ((Candidate1 - Candidate2)[1] < tol && (Candidate1 - Candidate2)[1] > -tol) &&
        ((Candidate1 - Candidate2)[2] < tol && (Candidate1 - Candidate2)[2] > -tol))
    {
        intersection = Candidate1;
        return true;
    }
    else
        return false;
}

bool isPointInPolygon(const Eigen::Vector3d& point,
                      const Eigen::MatrixXd& Fracture)
{
    double tol = 1e-10;
    int numVertices = Fracture.cols();

    if (numVertices < 4)
    {
        return false;
    }

    Eigen::Vector3d AB = Fracture.col(1) - Fracture.col(0);
    Eigen::Vector3d AD = Fracture.col(3) - Fracture.col(0); //ricorda che hai cambiato il tre in due
    Eigen::Vector3d AP = point - Fracture.col(0);

    double dotAB_AB = AB.dot(AB);
    double dotAD_AD = AD.dot(AD);
    double dotAP_AB = AP.dot(AB);
    double dotAP_AD = AP.dot(AD);

    double lambda1 = dotAP_AB / dotAB_AB;
    double lambda2 = dotAP_AD / dotAD_AD;

    return (lambda1 >= -tol && lambda1 <= 1.0+tol && lambda2 >= -tol && lambda2 <= 1+tol);
}

bool isPointInPolygonOK(const Eigen::Vector3d& point,
                        const Data::Fract& Fracture)
{
    bool PointInPolygon = true;
    for(unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        Eigen::Vector3d vector1 = Fracture.vertices.col((i+1) % Fracture.vertices.cols()) - Fracture.vertices.col(i);
        Eigen::Vector3d vector2 = point - Fracture.vertices.col(i);
        Eigen::Vector3d CrossProduct = vector1.cross(vector2);
        if(CrossProduct.dot(Fracture.normals) < -tol )
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
        if(isPointInPolygon(FirstFracture.vertices.col(i), SecondFracture.vertices))
            candidatePoints.push_back(FirstFracture.vertices.col(i));
    }
    for (unsigned int i = 0; i < SecondFracture.vertices.cols(); i++)
    {
        if(isPointInPolygon(SecondFracture.vertices.col(i), FirstFracture.vertices))
            candidatePoints.push_back(SecondFracture.vertices.col(i));
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
    //chiedi  se escludere caso in cui l'intersezione è un solo punto
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
    //Eigen::Vector3d ExtremeTrace;
    FractureOperations::findPosition(FirstFracture, t, P, CandidatePoints);
    FractureOperations::findPosition(SecondFracture, t, P, CandidatePoints);

    std::vector<Eigen::VectorXd> extremePoints;
    for (unsigned int i = 0; i < CandidatePoints.size(); i++)
    {
        if (isPointInPolygon(CandidatePoints[i], FirstFracture.vertices) && isPointInPolygon(CandidatePoints[i], SecondFracture.vertices))
            extremePoints.push_back(CandidatePoints[i]);
    }

    if (extremePoints.size() != 0)
    {
        std::vector<Eigen::VectorXd> extremePointsOK;
        extremePointsOK.push_back(extremePoints[0]);

        for (unsigned int i = 1; i < extremePoints.size(); i++)
        {
            if (extremePoints[i] != extremePoints[i-1])
                extremePointsOK.push_back(extremePoints[i]);
        }

        foundTrace.FractureIds[0] = FirstFracture.FractId;
        foundTrace.FractureIds[1] = SecondFracture.FractId;

        foundTrace.ExtremesCoord[0] = extremePointsOK[0];
        foundTrace.ExtremesCoord[1] = extremePointsOK[1];

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
        { //cambio <= con >= per ordinare in modo decrescente e invece di prendere solo l'elemento ne guardo la lunghezza
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
    { //va aggiunto il controllo perché se data.size()=0 avrei come right -1 (ma size_t non può essere negativo)
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

    if ( ((Candidate1 - Candidate2)[0] < tol && (Candidate1 - Candidate2)[0] > -tol) &&
        ((Candidate1 - Candidate2)[1] < tol && (Candidate1 - Candidate2)[1] > -tol) &&
        ((Candidate1 - Candidate2)[2] < tol && (Candidate1 - Candidate2)[2] > -tol))
    {
        Solution = Candidate1;
        Flag = true;
    }

    return Flag;
}

/*
bool ConsideringTollerance(const Data::Fract& CurrentPolygon,
                           Eigen::MatrixXd& subpol1,
                           Eigen::MatrixXd& subpol2,
                           std::queue<Data::Fract>& AllSubPolygons)
{
    for (unsigned int i = 0; i < subpol1.cols(); i++)
        if((subpol1.col(i) - subpol1.col((i - 1) % subpol1.cols())).norm() < tol)
            return false;

    for (unsigned int i = 0; i < subpol2.cols(); i++)
        if((subpol2.col(i) - subpol2.col((i - 1) % subpol2.cols())).norm() < tol)
            return false;

    Data::Fract subpolygon;
    subpolygon.vertices = subpol1;
    subpolygon.passingTracesId = CurrentPolygon.passingTracesId;
    //subpolygon.notPassingTracesId= Fracture.notPassingTracesId;
    subpolygon.normals = CurrentPolygon.normals;
    //aggiungo alla fine della lista i sottopoligoni da analizzare
    //std::cout << SubPolygon << std::endl;
    //std::cout << std::endl;
    AllSubPolygons.push(subpolygon);

    Data::Fract secondsubpolygon;
    secondsubpolygon.vertices = subpol2;
    secondsubpolygon.passingTracesId = CurrentPolygon.passingTracesId;
    //subpolygon.notPassingTracesId= Fracture.notPassingTracesId;
    secondsubpolygon.normals = CurrentPolygon.normals;
    //aggiungo alla fine della lista i sottopoligoni da analizzare
    //std::cout << SubPolygon << std::endl;
    //std::cout << std::endl;
    AllSubPolygons.push(secondsubpolygon);

    return true;
}
*/

bool MakeCuts(std::list<unsigned int>& AllTraces,
              std::vector<Data::Trace>& traces,
              PolygonalMeshLibrary::PolygonalMesh& PolygonMesh,
              std::queue<Data::Fract>& AllSubPolygons)
{
    Data::Fract CurrentPolygon;

    if (AllTraces.size() == 0)
    {
        std::cout << "bisogna salvare il sottopoligono " << std::endl;

        while(!AllSubPolygons.empty())
        {
            PolygonalMeshLibrary::CreateMesh(AllSubPolygons.front(), PolygonMesh);
            AllSubPolygons.pop();
        }
        return true;
    }
    else
    {
        if(AllSubPolygons.size() == 0 )
        {
            std::cerr << "ERROR: AllSubPolygons.size() == 0 " << std::endl;
            return false;
        }
        else
        {
            CurrentPolygon = AllSubPolygons.front();
        }

        //devo selezionare la traccia più lunga fra quelle del sottopoligono
        Data::Trace CurrentTrace;
        //se metto reference vicino ad auto mi fa errore
        bool FineRicorsione = true;
        for (auto it = AllTraces.begin(); it != AllTraces.end(); it++)
        {
            if(FineRicorsione)
            {
                CurrentTrace = traces[*it];
                //std::cout << "TraceId " << CurrentTrace.TraceId << std::endl;
                // vedere come gestire casi in cui la traccia viene divisa fra due sotto poligoni e chiedere se bisogna farlo
                if (FractureOperations::isPointInPolygonOK(CurrentTrace.ExtremesCoord[0], CurrentPolygon) ||
                    FractureOperations::isPointInPolygonOK(CurrentTrace.ExtremesCoord[1], CurrentPolygon))
                {
                        /*quando la traccia è interna al sottopoligono
                     * FineRicorsione diventa false ed esco dal ciclo
                     * altrimenti controllo se la traccia successiva è interns al poligonp
                     * Se non trovo questa traccia passo al poligonp suvccessivo
                     * FineRicorsione diventa true
                    */
                    FineRicorsione = false;
                }
            }
        }

        std::cout << CurrentTrace.TraceId << " " << FineRicorsione << std::endl;

        if(FineRicorsione == true)
        {
            // bisogna salvare il sotto poligono e uscire dalla funzione
            std::cout << " bisogna salvare il sotto poligono e passare al prossimo sootopoligono" << std::endl;
            PolygonalMeshLibrary::CreateMesh(CurrentPolygon,
                                             PolygonMesh);
            // Rimuovi il primo sottopoligono analizzato dalla lista
            AllSubPolygons.pop();

            //richiamo makecut per ogni sotto pologono
            PolygonalMeshLibrary::MakeCuts(AllTraces,
                                           traces,
                                           PolygonMesh,
                                           AllSubPolygons);
        }

        //estraggo gli estremi e individuo la direzione sulla quale giace la traccia
        Eigen::Vector3d FirstExtreme = CurrentTrace.ExtremesCoord[0];
        Eigen::Vector3d SecondExtreme = CurrentTrace.ExtremesCoord[1];
        Eigen::Vector3d Direction = SecondExtreme - FirstExtreme;
        //per ogni vertice del poligono controllo la sua posizione reciproca rispetto alla traccia
        //inizio a creare un std::vector ma fose meglio array(bisogna controllare costi computazionali) sappiamo dimensioni a priori??

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


            if(FractureOperations::isPointOnEdge(FirstExtreme, CurrentPolygon.vertices.col(i), CurrentPolygon.vertices.col((i - 1) % CurrentPolygon.vertices.cols())))
            {
                estremiTracce.push_back(FirstExtreme);
                Eigen::Vector3d congiungente = CurrentPolygon.vertices.col(i) - FirstExtreme;
                Eigen::Vector3d v= Direction.cross(congiungente);
                if (v.dot(CurrentPolygon.normals)>-tol) //controlla > tol oppure >= -tol
                    CurrentCheck = true;
            }
            else if(FractureOperations::isPointOnEdge(SecondExtreme, CurrentPolygon.vertices.col(i), CurrentPolygon.vertices.col((i - 1) % CurrentPolygon.vertices.cols())))
            {
                estremiTracce.push_back(SecondExtreme);
                Eigen::Vector3d congiungente = CurrentPolygon.vertices.col(i) - FirstExtreme;
                Eigen::Vector3d v= Direction.cross(congiungente);
                if (v.dot(CurrentPolygon.normals)>-tol) //controlla > tol oppure >= -tol
                    CurrentCheck = true;
            }
            else
            {
                Eigen::Vector3d congiungente = CurrentPolygon.vertices.col(i) - FirstExtreme;
                Eigen::Vector3d v= Direction.cross(congiungente);
                if (v.dot(CurrentPolygon.normals)>-tol) //controlla > tol oppure >= -tol
                    CurrentCheck = true;
                if(PreviousCheck != CurrentCheck)
                {
                    //solve sistem
                    Eigen::Vector3d Solution;
                    if(PolygonalMeshLibrary::SolveSystem(Direction,
                                                          CurrentPolygon.vertices.col(i),
                                                          CurrentPolygon.vertices.col((i - 1) % CurrentPolygon.vertices.cols()),
                                                          FirstExtreme,
                                                          Solution))
                    {
                        estremiTracce.push_back(Solution);
                    }
                }
            }
            PreviousCheck = CurrentCheck;
        }

        for (auto& elem : estremiTracce)
            std::cout <<elem.transpose() << std::endl;

        if((estremiTracce.front() - estremiTracce.back()).norm() < tol )
        {
            std::cout << " estremi tracce è un punto " << std::endl;
            //AllSubPolygons.pop();
            /*
            PolygonalMeshLibrary::MakeCuts(AllTraces,
                                           traces,
                                           PolygonalMesh,
                                           AllSubPolygons);
    */
        }
        FirstExtreme = estremiTracce.front();
        Direction = estremiTracce.back() - FirstExtreme;
        std::vector<Eigen::Vector3d> FirstSide;

        unsigned int i = 0;
        congiungente = CurrentPolygon.vertices.col(i) - FirstExtreme;

        v= Direction.cross(congiungente);

        while (v.dot(CurrentPolygon.normals)> tol )
        {
            FirstSide.push_back(CurrentPolygon.vertices.col(i++));
            congiungente = CurrentPolygon.vertices.col(i) - FirstExtreme;
            v= Direction.cross(congiungente);
        }
        //forse prima back e poi front
        FirstSide.push_back(estremiTracce.front());
        FirstSide.push_back(estremiTracce.back());

        for(unsigned int j = i; j < CurrentPolygon.vertices.cols(); j++)
        {
            Eigen::Vector3d congiungente = CurrentPolygon.vertices.col(j) - FirstExtreme;
            Eigen::Vector3d v= Direction.cross(congiungente);
            if (v.dot(CurrentPolygon.normals)>tol)
            {
                FirstSide.push_back(CurrentPolygon.vertices.col(j));
                //se la normale è parallela sta da un lato
                //FirstSide.push_back(CurrentPolygon.vertices.col(i));
                //std::cout << CurrentPolygon.vertices.col(i).transpose() << std::endl;
            }
        }

        std::vector<Eigen::Vector3d> SecondSide;
        congiungente = CurrentPolygon.vertices.col(0) - FirstExtreme;
        v= Direction.cross(congiungente);

        unsigned int k = 0;
        while (v.dot(CurrentPolygon.normals) < - tol )
        {
            SecondSide.push_back(CurrentPolygon.vertices.col(k++));
            congiungente = CurrentPolygon.vertices.col(k) - FirstExtreme;
            v= Direction.cross(congiungente);
        }
        //forse prima back e poi front
        SecondSide.push_back(estremiTracce.back());
        SecondSide.push_back(estremiTracce.front());

        for(unsigned int j =k; j < CurrentPolygon.vertices.cols(); j++)
        {
            Eigen::Vector3d congiungente = CurrentPolygon.vertices.col(j) - FirstExtreme;
            Eigen::Vector3d v= Direction.cross(congiungente);
            if (v.dot(CurrentPolygon.normals)< - tol)
            {
                SecondSide.push_back(CurrentPolygon.vertices.col(j));
                //se la normale è parallela sta da un lato
                //FirstSide.push_back(CurrentPolygon.vertices.col(i));
                //std::cout << CurrentPolygon.vertices.col(i).transpose() << std::endl;
            }
        }

        //se uno dei due poligoni è degenere, vuol dire la traccia è su un lato => l'altro poligono coincide con
        //current polygon pertanto non lo salvo
        if(FirstSide.size() > 2 && SecondSide.size() >2)
        {
            //copio std::vector in Eigen::Matrix
            Eigen::MatrixXd SubPolygon(3,FirstSide.size());

            for(unsigned int k = 0; k < FirstSide.size(); k++)
            {
                SubPolygon.col(k) = FirstSide[k];
            }

            Data::Fract subpolygon;
            subpolygon.vertices = SubPolygon;
            subpolygon.passingTracesId = CurrentPolygon.passingTracesId;
            //subpolygon.notPassingTracesId= Fracture.notPassingTracesId;
            subpolygon.normals = CurrentPolygon.normals;
            //aggiungo alla fine della lista i sottopoligoni da analizzare
            //std::cout << SubPolygon << std::endl;
            //std::cout << std::endl;
            AllSubPolygons.push(subpolygon);

            //copio std::vector in Eigen::Matrix
            Eigen::MatrixXd SecondSubPolygon(3,SecondSide.size());
            //std::cout << "oke" << std::endl;

            for(unsigned int k = 0; k < SecondSide.size(); k++)
            {
                SecondSubPolygon.col(k) = SecondSide[k];
            }

            Data::Fract secondsubpolygon;
            secondsubpolygon.vertices = SecondSubPolygon;
            secondsubpolygon.passingTracesId = CurrentPolygon.passingTracesId;
            //subpolygon.notPassingTracesId= Fracture.notPassingTracesId;
            secondsubpolygon.normals = CurrentPolygon.normals;
            //aggiungo alla fine della lista i sottopoligoni da analizzare
            //std::cout << SubPolygon << std::endl;
            //std::cout << std::endl;
            AllSubPolygons.push(secondsubpolygon);

            /*
            if (!PolygonalMeshLibrary::ConsideringTollerance(CurrentPolygon,
                                                             SubPolygon,
                                                             SecondSubPolygon,
                                                             AllSubPolygons))
            {
                std::cout << "non ce la faccio più" << std::endl;
            }
            */
        }
        else
        {
            //sicuramnete è sbagliato
            AllTraces.remove(CurrentTrace.TraceId);
            PolygonalMeshLibrary::MakeCuts(AllTraces,
                                           traces,
                                           PolygonMesh,
                                           AllSubPolygons);
        }

        //std::cout << SecondSide.size() << std::endl;

        //elimino la traccia considerata:
        if (FractureOperations::isPointInPolygonOK(CurrentTrace.ExtremesCoord[0], CurrentPolygon) &&
            FractureOperations::isPointInPolygonOK(CurrentTrace.ExtremesCoord[1], CurrentPolygon))
        {
            AllTraces.remove(CurrentTrace.TraceId);
        }
        else if (FractureOperations::isPointInPolygonOK(CurrentTrace.ExtremesCoord[0], CurrentPolygon))
        {

            std::cout << CurrentPolygon.vertices << std::endl;
            std::cout << std::endl;

            std::cout << CurrentTrace.ExtremesCoord[0].transpose() << std::endl;
            std::cout << CurrentTrace.ExtremesCoord[1].transpose() << std::endl;
            std::cout << std::endl;

            std::cout << estremiTracce.back().transpose() << std::endl;
            std::cout << estremiTracce.front().transpose() << std::endl;
            std::cout << std::endl;


            //vedi bene front e back, perchè secondo me è al contrario
            //aggiungi tol
            if((CurrentTrace.ExtremesCoord[1] - estremiTracce.back()).norm() < (CurrentTrace.ExtremesCoord[1] - estremiTracce.front()).norm() )
                traces[CurrentTrace.TraceId].ExtremesCoord[0] = estremiTracce.back();
            else
                traces[CurrentTrace.TraceId].ExtremesCoord[0] = estremiTracce.front();

            if ((CurrentTrace.ExtremesCoord[0] - CurrentTrace.ExtremesCoord[1]).norm() < tol)
            {
                AllTraces.remove(CurrentTrace.TraceId);
            }

            std::cout << traces[CurrentTrace.TraceId].ExtremesCoord[0].transpose() << std::endl;
            std::cout << traces[CurrentTrace.TraceId].ExtremesCoord[1].transpose() << std::endl;
            //std::cout << std::endl;
        }
        else if (FractureOperations::isPointInPolygonOK(CurrentTrace.ExtremesCoord[1], CurrentPolygon))
        {
            if((CurrentTrace.ExtremesCoord[0] - estremiTracce.back()).norm() < (CurrentTrace.ExtremesCoord[0] - estremiTracce.front()).norm() )
                traces[CurrentTrace.TraceId].ExtremesCoord[1] = estremiTracce.back();
            else
                traces[CurrentTrace.TraceId].ExtremesCoord[1] = estremiTracce.front();

            if ((CurrentTrace.ExtremesCoord[0] - CurrentTrace.ExtremesCoord[1]).norm() < tol)
            {
                AllTraces.remove(CurrentTrace.TraceId);

            }
            //std::cout << "problem qua" << std::endl;
        }


        // Rimuovi il primo sottopoligono analizzato dalla lista
        AllSubPolygons.pop();

        //richiamo makecut per ogni sotto pologono
        //la chiamo una sola volta perchè uso sempre il primo sottopoligono

        PolygonalMeshLibrary::MakeCuts(AllTraces,
                                       traces,
                                       PolygonMesh,
                                       AllSubPolygons);
    }

    return true;
}

void CreateMesh(const Data::Fract& Fracture,
                PolygonalMeshLibrary::PolygonalMesh& PolygonMesh)
{
    std::cout << Fracture.vertices << std::endl;
    std::cout << std::endl;

    std::vector<unsigned int> IdVectices;
    //salvo i vertici
    for(unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        auto variabile = std::find(PolygonMesh.coord0DsCell.begin(), PolygonMesh.coord0DsCell.end(), Fracture.vertices.col(i));
        if(variabile == PolygonMesh.coord0DsCell.end())
        {
            //se punta a PolygonalMesh.coord0DsCell.end() cevo aggiungere
            PolygonMesh.coord0DsCell.push_back(Fracture.vertices.col(i));
            IdVectices.push_back(PolygonMesh.coord0DsCell.size()-1);
        }
        else
        {
            unsigned int position = std::distance(PolygonMesh.coord0DsCell.begin(), variabile);
            IdVectices.push_back(position);
        }
    }


    /* ho pensato una cvosa ,a non so se può avere senso (fammi sap)
     * noi per ogni polignono, cicliamo sui vertici e  salviamo l'id in Polygonal mesh
     * poi per ogni polignono, cicliamo sui lati e  salviamo le coordinate in Polygonal mesh (dovremmo salvare gli Id)
     * se facessimo così:
     * per ogni polignono, cicliamo sui vertici e  salviamo l'id in Polygonal mesh e li salviamo in un vector o lista
     * tramite l'indice di PolygonalMesh.coord0DsCell
     * poi cicliamo for(unsigned int i = 0; i < Fracture.vertices.cols(); i++) e con i accediamo all'indice
     */
     // secondo me, l'ideale sarebbe usare un std::map a cui accedere con chiave l'id perchè è più basse il costo del find ma non ho capito come si fa a cercare per valore


    std::vector<unsigned int> IdEdge;

    for(unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        std::array<unsigned int, 2> edge1 = {IdVectices[i], IdVectices[(i-1) % Fracture.vertices.cols()]};
        std::array<unsigned int, 2> edge2 = {IdVectices[(i-1) % Fracture.vertices.cols()], IdVectices[i]};
        auto variabile1 = std::find(PolygonMesh.coord1DsCell_Id0DS.begin(), PolygonMesh.coord1DsCell_Id0DS.end(), edge1);
        auto variabile2 = std::find(PolygonMesh.coord1DsCell_Id0DS.begin(), PolygonMesh.coord1DsCell_Id0DS.end(), edge2);
        if(variabile1 == PolygonMesh.coord1DsCell_Id0DS.end() && variabile2 == PolygonMesh.coord1DsCell_Id0DS.end() )
        {
            //se punta a PolygonalMesh.coord0DsCell.end() cevo aggiungere
            PolygonMesh.coord1DsCell_Id0DS.push_back(edge1);
            IdEdge.push_back(PolygonMesh.coord0DsCell.size());
        }
        else
        {
            if (variabile1 == PolygonMesh.coord1DsCell_Id0DS.end())
            {
            unsigned int position = std::distance(PolygonMesh.coord1DsCell_Id0DS.begin(), variabile2);
            IdEdge.push_back(position-1);
            }
            else
            {
                unsigned int position = std::distance(PolygonMesh.coord1DsCell_Id0DS.begin(), variabile1);
                IdEdge.push_back(position-1);
            }
        }
    }

    PolygonMesh.Cell2DsVertices.push_back(IdVectices);
    PolygonMesh.Cell2DsEdges.push_back(IdEdge);
}
}




