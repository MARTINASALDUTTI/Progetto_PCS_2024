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
    //double tol = 1e-10;
 //   if (Fracture.cols() <3)
 //   {
  //      std::cerr << "non si tratta di un poligono" << std::endl;
  //  }
    //VEDI TRE: POSSIBILE ERRORE
    if(Fracture.cols() > 3 )
    {
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
    return false;
}

bool isPointInPolygonCorrectversion (const Eigen::Vector3d& point,
                                    const Eigen::MatrixXd& Fracture)
{
    double tol = 1e-10;
    int numVertices = Fracture.cols();

    if (numVertices < 3)
    {
        std::cerr << "Il poligono deve avere almeno 3 vertici." << std::endl;
        return false;
    }

    Eigen::Vector3d p0 = Fracture.col(0);
    for (int i = 1; i <= numVertices; ++i)
    {
        Eigen::Vector3d p1 = Fracture.col(i % numVertices);
        Eigen::Vector3d p2 = Fracture.col((i + 1) % numVertices);

        Eigen::Vector3d edge1 = p1 - p0;
        Eigen::Vector3d edge2 = p2 - p0;
        Eigen::Vector3d edge3 = p2 - p1;

        Eigen::Vector3d normal = edge1.cross(edge2).normalized();

        Eigen::Vector3d AP = point - p0;
        Eigen::Vector3d BP = point - p1;

        double dotAP_normal = AP.dot(normal);
        double dotBP_normal = BP.dot(normal);

        if (dotAP_normal < -tol || dotBP_normal < -tol)
            return false;

        p0 = p1;
    }

    return true;
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
Eigen::MatrixXd SavingPoint(const std::vector<Eigen::Vector3d>& PointsList,
                            const Eigen::Vector3d& Normal, Eigen::MatrixXd &MatrixPoint)
{
    Eigen::MatrixXd MatrixPoint(3, PointsList.size());

    for(unsigned int k = 0; k < PointsList.size(); k++)
    {
        //....
    }
    return MatrixPoint;
}
*/

bool MakeCuts(std::list<unsigned int>& AllTraces,
              const std::vector<Data::Trace>& traces,
              PolygonalMeshLibrary::PolygonalMesh& PolygonalMesh,
              std::queue<Data::Fract>& AllSubPolygons)
{
    Eigen::Vector3d Candidate1;
    Eigen::Vector3d Candidate2;
    Data::Fract Fracture;

    if (AllTraces.size() == 0)
    {
        //std::cout << "bisogna salvare il sottopoligono " << std::endl;
         // non so se bisogna aggiungere ciclare per tutti (forse sì)
        //fine della ricorsione
        while(!AllSubPolygons.empty())
        {
            PolygonalMeshLibrary::CreateMesh(AllSubPolygons.front(), PolygonalMesh);
            AllSubPolygons.pop();
        }
        //Fracture = AllSubPolygons.front();
        //PolygonalMeshLibrary::CreateMesh(Fracture,PolygonalMesh);

        return false;
    }
    else
    {
        if (AllSubPolygons.size() == 0)
        {
            std::cout << " AllTraces.size() " << AllTraces.size() << std::endl;
            return false;
        }
        else
            Fracture = AllSubPolygons.front();

        //seleziono la traccia più lunga
        //forse meglio lista
        Data::Trace CurrentTrace;
        //se metto reference vicino ad auto mi fa errore
        bool FineRicorsione = true;
        for (auto it = AllTraces.begin(); it != AllTraces.end(); ++it)
        {
            CurrentTrace = traces[*it];
            //std::cout << "TraceId " << CurrentTrace.TraceId << std::endl;
            // vedere come gestire casi in cui la traccia viene divisa fra due sotto poligoni e chiedere se bisogna farlo

            if (FractureOperations::isPointInPolygonCorrectversion(CurrentTrace.ExtremesCoord[0], Fracture.vertices))
            {
                /*quando la traccia è interna al sottopoligono
                 * FineRicorsione diventa false ed esco dal ciclo
                 * altrimenti controllo se la traccia successiva è interns al poligonp
                 * Se non trovo questa traccia passo al poligonp suvccessivo
                 * FineRicorsione diventa true
                */
                FineRicorsione = false;
                //vedere se usare break (pericoloso o un if su fine ricorsione)
                break;
            }
            //capire cosa fare se esce un segmento e perche
        }

        //std::cout << "FineRicorsione " << FineRicorsione << std::endl;
        if(FineRicorsione == true)
        {
            // bisogna salvare il sotto poligono e uscire dalla funzione
            std::cout << " bisogna salvare il sotto poligono e passare al prossimo sootopoligono" << std::endl;
            PolygonalMeshLibrary::CreateMesh(Fracture,PolygonalMesh);
            // Rimuovi il primo sottopoligono analizzato dalla lista
            AllSubPolygons.pop();

            //richiamo makecut per ogni sotto pologono
            PolygonalMeshLibrary::MakeCuts(AllTraces,
                                           traces,
                                           PolygonalMesh,
                                           AllSubPolygons);
            return false;
        }
        //estraggo gli estremi e individuo la direzione sulla quale giace la traccia
        Eigen::Vector3d FirstExtreme = CurrentTrace.ExtremesCoord[0];
        Eigen::Vector3d SecondExtreme = CurrentTrace.ExtremesCoord[1];
        Eigen::Vector3d Direction = SecondExtreme - FirstExtreme;
        //per ogni vertice del poligono controllo la sua posizione reciproca rispetto alla traccia
        //inizio a creare un std::vector ma fose meglio array(bisogna controllare costi computazionali) sappiamo dimensioni a priori??
        std::vector<Eigen::Vector3d> FirstSide;
        std::vector<Eigen::Vector3d> SecondSide;
        //se la traccia è passante già ho gli estremi che appartegono ad entrambi i poligoni
        bool Passing = false;
        //tramite tip
        auto variabile = std::find(Fracture.passingTracesId.begin(), Fracture.passingTracesId.end(), CurrentTrace.TraceId);
        if(variabile != Fracture.passingTracesId.end())
        {
            //creo un bool per indicare che la traccia è passante
            Passing = true;
        }

        //std::cout << "Passing " << Passing << std::endl;
        //se non è passante devo prolungare
        if(Passing == false)
        {

            std::vector<Eigen::Vector3d> estremiTracce;

            bool PreviousCheck = false;
            Eigen::Vector3d congiungente = Fracture.vertices.col(Fracture.vertices.cols()-1) - FirstExtreme;
            Eigen::Vector3d v= Direction.cross(congiungente);
            if (v.dot(Fracture.normals)>-tol)
                PreviousCheck = true;
            for (unsigned int i = 0; i < Fracture.vertices.cols(); i++)
            {
                /*per ogni lato ho tre opzioni
                 * 1 il primo estremo è sul lato
                 * 2 il secondo estremo è sul lato
                 * 3 nessuno dei due estremi è sul lato, controllo se devo prolungare
                 */
                bool CurrentCheck = false;

                if(FractureOperations::isPointOnEdge(FirstExtreme, Fracture.vertices.col(i), Fracture.vertices.col((i - 1) % Fracture.vertices.cols())))
                {
                    estremiTracce.push_back(FirstExtreme);
                    Eigen::Vector3d congiungente = Fracture.vertices.col(i) - FirstExtreme;
                    Eigen::Vector3d v= Direction.cross(congiungente);
                    if (v.dot(Fracture.normals)>-tol) //controlla > tol oppure >= -tol
                        CurrentCheck = true;
                }
                else if(FractureOperations::isPointOnEdge(SecondExtreme, Fracture.vertices.col(i), Fracture.vertices.col((i - 1) % Fracture.vertices.cols())))
                {
                    estremiTracce.push_back(SecondExtreme);
                    Eigen::Vector3d congiungente = Fracture.vertices.col(i) - FirstExtreme;
                    Eigen::Vector3d v= Direction.cross(congiungente);
                    if (v.dot(Fracture.normals)>-tol) //controlla > tol oppure >= -tol
                        CurrentCheck = true;
                }
                else
                {
                    Eigen::Vector3d congiungente = Fracture.vertices.col(i) - FirstExtreme;
                    Eigen::Vector3d v= Direction.cross(congiungente);
                    if (v.dot(Fracture.normals)>-tol) //controlla > tol oppure >= -tol
                        CurrentCheck = true;
                    if(PreviousCheck != CurrentCheck)
                    {
                        //solve sistem
                        Eigen::Vector3d Solution;
                        if(PolygonalMeshLibrary::SolveSystem(Direction,
                                                             Fracture.vertices.col(i),
                                                             Fracture.vertices.col((i - 1) % Fracture.vertices.cols()),
                                                             FirstExtreme,
                                                             Solution))
                        {
                            estremiTracce.push_back(Solution);
                        }
                    }
                    else if(((Direction.dot(Fracture.vertices.col(i) - FirstExtreme) >= -tol) || (Direction.dot(Fracture.vertices.col(i) - FirstExtreme) <= tol)) &&
                            ((Direction.dot(Fracture.vertices.col((i - 1) % Fracture.vertices.cols()) - FirstExtreme) >= -tol) ||
                             (Direction.dot(Fracture.vertices.col((i - 1) % Fracture.vertices.cols()) - FirstExtreme) <= tol)))
                    {
                        //solve sistem
                        Eigen::Vector3d Solution;
                        if(PolygonalMeshLibrary::SolveSystem(Direction,
                                                             Fracture.vertices.col(i),
                                                             Fracture.vertices.col((i - 1) % Fracture.vertices.cols()),
                                                             FirstExtreme,
                                                             Solution))
                        {
                            estremiTracce.push_back(Solution);
                        }
                    }
                }
                PreviousCheck = CurrentCheck;
            }
            //aggiorno FirstExtreme e SecondExtreme con i prolungamenti
            //estremi tracce ha necessariamente due elementi
            FirstExtreme = estremiTracce[0];
            SecondExtreme = estremiTracce[1];
        }
         //ho prolungato quindi nel caso di tracce passanti ho gli estremi aaltrimenti i prolungamenti
        FirstSide.push_back(FirstExtreme);
        SecondSide.push_back(FirstExtreme);

        bool AnotherCheck = false;
        for (unsigned int i = 0; i < Fracture.vertices.cols(); i++)
        {
            Eigen::Vector3d congiungente = Fracture.vertices.col(i) - FirstExtreme;
            Eigen::Vector3d v= Direction.cross(congiungente);
            if (v.dot(Fracture.normals)>tol)
            {
                //se la normale è parallela sta da un lato
                FirstSide.push_back(Fracture.vertices.col(i));
            }
            else if (v.dot(Fracture.normals)<-tol)
            {
                SecondSide.push_back(Fracture.vertices.col(i));
            }
            else if (v.dot(Fracture.normals)>-tol && v.dot(Fracture.normals)<tol)
            {
                //caso limite: il prodotto scalare fa 0
                std::cout << "caso limite" << std::endl;
                AnotherCheck = true;
                break;
            }
        }

        //l'altro controllo: gestione caso limite ( valido solo per quadrati)
        if(AnotherCheck == true)
        {
            Eigen::Vector3d edge = Fracture.vertices.col(0) - Fracture.vertices.col(1);
            if(Direction.cross(edge).norm()> -tol && Direction.cross(edge).norm() < tol)
            {
                //std::cout << "caso pippo 1" << std::endl;

                FirstSide.push_back(Fracture.vertices.col(0));
                FirstSide.push_back(Fracture.vertices.col(1));
                SecondSide.push_back(Fracture.vertices.col(2));
                SecondSide.push_back(Fracture.vertices.col(3));
            }
            else
            {
                //std::cout << "caso pippo 2" << std::endl;

                FirstSide.push_back(Fracture.vertices.col(1));
                FirstSide.push_back(Fracture.vertices.col(2));
                SecondSide.push_back(Fracture.vertices.col(3));
                SecondSide.push_back(Fracture.vertices.col(0));
            }
        }

        //per rispettare il verso antiorario
        FirstSide.push_back(SecondExtreme);
        SecondSide.push_back(SecondExtreme);

        //bisogna capire perchè escono poligoni da 2
        if(FirstSide.size() > 2)
        {
            //Eigen::MatrixXd FirstSubPolygon = SavingPoint(FirstSide, Fracture.normals);
            //copio std::vector in Eigen::Matrix
            Eigen::MatrixXd FirstSubPolygon(3,FirstSide.size());

            for(unsigned int k = 0; k < FirstSide.size(); k++)
            {
                FirstSubPolygon.col(k) = FirstSide[k];
            }
            //std::cout << FirstSubPolygon << std::endl;
            Data::Fract subpolygonuno;
            subpolygonuno.vertices = FirstSubPolygon;
            subpolygonuno.passingTracesId = Fracture.passingTracesId;
            subpolygonuno.notPassingTracesId= Fracture.notPassingTracesId;
            subpolygonuno.normals = Fracture.normals;
            //aggiungo alla fine della lista i sottopoligoni da analizzare

            AllSubPolygons.push(subpolygonuno);
        }

        if(SecondSide.size() > 2)
        {
            Eigen::MatrixXd SecondSubPolygon(3,SecondSide.size());

            for(unsigned int k = 0; k < SecondSide.size(); k++)
            {
                SecondSubPolygon.col(k) = SecondSide[k];
            }
            //std::cout << SecondSubPolygon << std::endl;
            Data::Fract subpolygondue;
            subpolygondue.vertices = SecondSubPolygon;
            subpolygondue.passingTracesId = Fracture.passingTracesId;
            subpolygondue.notPassingTracesId= Fracture.notPassingTracesId;
            subpolygondue.normals = Fracture.normals;
            //aggiungo alla fine della lista i sottopoligoni da analizzare

            AllSubPolygons.push(subpolygondue);
        }

        /*
        std::cout << " AllSubPolygons " << std::endl;
        for (auto it = AllSubPolygons.begin(); it != AllSubPolygons.end(); ++it)
        {
            Fracture = *it;
            PolygonalMeshLibrary::CreateMesh(Fracture, PolygonalMesh);
        }
        */
        //elimino la traccia considerata
        AllTraces.remove(CurrentTrace.TraceId); //remove elimina tutte le occorrenze... NO PROBLEM: nell nostra lista non ci sono elementi ripetuti

        // Rimuovi il primo sottopoligono analizzato dalla lista
        AllSubPolygons.pop();

        //richiamo makecut per ogni sotto pologono
        //la chiamo una sola volta perchè uso sempre il primo sottopoligono

        PolygonalMeshLibrary::MakeCuts(AllTraces,
                                       traces,
                                       PolygonalMesh,
                                       AllSubPolygons);
    }
    return true;
}

void CreateMesh(const Data::Fract& Fracture,
                PolygonalMeshLibrary::PolygonalMesh& PolygonalMesh)
{
    //std::cout << Fracture.vertices << std::endl;
    //std::cout << std::endl;

    std::vector<unsigned int> IdVectices;
    //salvo i vertici
    for(unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        auto variabile = std::find(PolygonalMesh.coord0DsCell.begin(), PolygonalMesh.coord0DsCell.end(), Fracture.vertices.col(i));
        if(variabile == PolygonalMesh.coord0DsCell.end())
        {
            //se punta a PolygonalMesh.coord0DsCell.end() cevo aggiungere
            PolygonalMesh.coord0DsCell.push_back(Fracture.vertices.col(i));
            IdVectices.push_back(PolygonalMesh.coord0DsCell.size()-1);
        }
        else
        {
            unsigned int position = std::distance(PolygonalMesh.coord0DsCell.begin(), variabile);
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
    /* con queto ciclo salvo le coordinate degli estremi ma in teoria non serve
    for(unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        std::array<Eigen::Vector3d, 2> edge1 = {Fracture.vertices.col(i), Fracture.vertices.col((i - 1) % Fracture.vertices.cols())};
        std::array<Eigen::Vector3d, 2> edge2 = {Fracture.vertices.col((i - 1) % Fracture.vertices.cols()), Fracture.vertices.col(i)};
        auto variabile1 = std::find(PolygonalMesh.coord1DsCell.begin(), PolygonalMesh.coord1DsCell.end(), edge1);
        auto variabile2 = std::find(PolygonalMesh.coord1DsCell.begin(), PolygonalMesh.coord1DsCell.end(), edge2);

        if(variabile1 == PolygonalMesh.coord1DsCell.end() && variabile2 == PolygonalMesh.coord1DsCell.end() )
        {
            //se punta a PolygonalMesh.coord0DsCell.end() cevo aggiungere
            PolygonalMesh.coord1DsCell.push_back(edge1);
        }
    }
    */

    std::vector<unsigned int> IdEdge;

    for(unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        std::array<unsigned int, 2> edge1 = {IdVectices[i], IdVectices[(i-1) % Fracture.vertices.cols()]};
        std::array<unsigned int, 2> edge2 = {IdVectices[(i-1) % Fracture.vertices.cols()], IdVectices[i]};
        auto variabile1 = std::find(PolygonalMesh.coord1DsCell_Id0DS.begin(), PolygonalMesh.coord1DsCell_Id0DS.end(), edge1);
        auto variabile2 = std::find(PolygonalMesh.coord1DsCell_Id0DS.begin(), PolygonalMesh.coord1DsCell_Id0DS.end(), edge2);
        if(variabile1 == PolygonalMesh.coord1DsCell_Id0DS.end() && variabile2 == PolygonalMesh.coord1DsCell_Id0DS.end() )
        {
            //se punta a PolygonalMesh.coord0DsCell.end() cevo aggiungere
            PolygonalMesh.coord1DsCell_Id0DS.push_back(edge1);
            IdEdge.push_back(PolygonalMesh.coord0DsCell.size());
        }
        else
        {
            if (variabile1 == PolygonalMesh.coord1DsCell_Id0DS.end())
            {
            unsigned int position = std::distance(PolygonalMesh.coord1DsCell_Id0DS.begin(), variabile2);
            IdEdge.push_back(position-1);
            }
            else
            {
                unsigned int position = std::distance(PolygonalMesh.coord1DsCell_Id0DS.begin(), variabile1);
                IdEdge.push_back(position-1);
            }
        }
    }

    PolygonalMesh.Cell2DsVertices.push_back(IdVectices);
    PolygonalMesh.Cell2DsEdges.push_back(IdEdge);

    /*
    for(unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {

        auto variabile1 = std::find(PolygonalMesh.coord1DsCell.begin(), PolygonalMesh.coord0DsCell.end(), Fracture.vertices.col(i));
        auto variabile2 = std::find(PolygonalMesh.coord1DsCell.begin(), PolygonalMesh.coord0DsCell.end(), Fracture.vertices.col(i));

        if(variabile == PolygonalMesh.coord0DsCell.end())
        {
            //se punta a PolygonalMesh.coord0DsCell.end() cevo aggiungere
            coord0DsCell.push_back(Fracture.vertices.col(i));
        }
    }

    std::unordered_map<Eigen::Vector3d*, unsigned int> coord0DsCellMap; //ha senso??
    //metodo esercitazione 5
    for(unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        const Eigen::Vector3d* vertice = &Fracture.vertices.col(i);

        //se la chiuave già esiste ret.second = false
        auto ret = PolygonalMesh.coord0DsCellMap.insert({vertice, {Cell0DId}});
        if(ret.second)
            Cell0DId++;
    }
    PolygonalMesh.coord0DsCellMap = coord0DsCellMap;
    */
    //una volta salvati i vertici, passo ai lati
    /*
    std::unordered_map<unsigned int, std::array<unsigned int, 2>> coord1DsCellMap;
    for(unsigned int i = 0; i < Fracture.vertices.cols(); i++)
    {
        unsigned int FirstId = coord0DsCellMap.find(Fracture.vertices.col(i)) ;
        unsigned int SecondId = coord0DsCellMap.find(Fracture.vertices.col(i) % Fracture.vertices.cols()) ;
        std::array<unsigned int, 2> value = {FirstId, SecondId};
        if (!containsValue(coord1DsCellMap, value))
        {
            coord1DsCellMap[(Cell1DId)++] = value; // Usa l'indice i come chiave
        }
    }
    PolygonalMesh.coord1DsCellMap = coord1DsCellMap;
    */
}
}




