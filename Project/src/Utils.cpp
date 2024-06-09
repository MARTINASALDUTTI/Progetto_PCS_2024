#include "fstream"
#include <iostream>
#include <sstream>
#include "Eigen/Eigen"
#include "iomanip"
#include <algorithm> // Per std::find

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

            if( Fractures[i].notPassingTracesId.size() != 0)
            {
                SortLibrary::Mergesort(Fractures[i].notPassingTracesId, Traces);
                for(auto& elem :Fractures[i].notPassingTracesId )
                {
                    outFile << elem << "; " << 0 << "; " << Traces[elem].length << "\n";
                }
            }
            if (Fractures[i].passingTracesId.size() != 0)
            {
                SortLibrary::Mergesort(Fractures[i].passingTracesId, Traces);
                for(auto& elem :Fractures[i].passingTracesId )
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

    Eigen::Vector3d AB = Fracture.col(1) - Fracture.col(0);
    Eigen::Vector3d AD = Fracture.col(3) - Fracture.col(0);
    Eigen::Vector3d AP = point - Fracture.col(0);

    double dotAB_AB = AB.dot(AB);
    double dotAD_AD = AD.dot(AD);
    double dotAP_AB = AP.dot(AB);
    double dotAP_AD = AP.dot(AD);

    double lambda1 = dotAP_AB / dotAB_AB;
    double lambda2 = dotAP_AD / dotAD_AD;

    return (lambda1 >= -tol && lambda1 <= 1.0+tol && lambda2 >= -tol && lambda2 <= 1+tol);
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

        foundTrace.Tips[0] = isTracePassing(FirstFracture.vertices, foundTrace.ExtremesCoord[0], foundTrace.ExtremesCoord[1]);
        foundTrace.Tips[1] = isTracePassing(SecondFracture.vertices, foundTrace.ExtremesCoord[0], foundTrace.ExtremesCoord[1]);

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

        foundTrace.Tips[0] = isTracePassing(FirstFracture.vertices, foundTrace.ExtremesCoord[0], foundTrace.ExtremesCoord[1]);
        foundTrace.Tips[1] = isTracePassing(SecondFracture.vertices, foundTrace.ExtremesCoord[0], foundTrace.ExtremesCoord[1]);

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
bool MakeCuts(const Data::Fract& Fracture,
              std::list<unsigned int>& AllTraces,
              const std::vector<Data::Trace>& traces,
              PolygonalMeshLibrary::PolygonalMesh& PolygonalMesh,
              unsigned int& Cell0DId,
              unsigned int& Cell1DId)
{
    Eigen::Vector3d Candidate1;
    Eigen::Vector3d Candidate2;

    if (AllTraces.size() == 0)
    {
        //fine della ricorsione
        std::cout << " bisogna salvare il sotto poligono" << std::endl;
        //PolygonalMeshLibrary::MakeMesh(Fracture,PolygonalMesh, Cell0DId, Cell1DId);
        return false;
    }
    else
    {
        //seleziono la traccia più lunga
        //forse meglio lista
        Data::Trace CurrentTrace;
        //se metto reference vicino ad auto mi fa errore
        for (auto it = AllTraces.begin(); it != AllTraces.end(); ++it)
        {
            CurrentTrace = traces[*it];
            if (FractureOperations::isPointInPolygon(CurrentTrace.ExtremesCoord[0], Fracture.vertices))
            {
                break;
            }
            else
            {
                // bisogna salvare il sotto poligono e uscire dalla funzione
                std::cout << " bisogna salvare il sotto poligono e uscire dalla funzione" << std::endl;
                //PolygonalMeshLibrary::MakeMesh(Fracture,PolygonalMesh, Cell0DId, Cell1DId);
                return false;
            }
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
        auto variabile = std::find(Fracture.passingTracesId.begin(), Fracture.passingTracesId.end(), CurrentTrace.TraceId);
        if(variabile != Fracture.passingTracesId.end())
        {
            //creo un bool per indicare che la traccia è passante
            Passing = true;
            std::cout << "passante " << std::endl;
        }
        else
            std::cout << " no passante " << std::endl;

        //se non è passante devo prolungare
        if(Passing == false)
        {
            std::vector<Eigen::Vector3d> estremiTracce;
            bool PreviousCheck = false;
            Eigen::Vector3d congiungente = Fracture.vertices.col(0) - FirstExtreme;
            Eigen::Vector3d v= Direction.cross(congiungente);
            if (v.dot(Fracture.normals)>tol)
                PreviousCheck = true;
            for (unsigned int i = 0; i < Fracture.vertices.cols(); i++)
            {
                /*per ogni lato ho tre opzioni
                 * 1 il primo estremo è sul lato
                 * 2 il secondo estremo è sul lato
                 * 3 nessuno dei due estremi è sul lato, controllo se devo prolungfare
                 */

                if(FractureOperations::isPointOnEdge(FirstExtreme, Fracture.vertices.col(i), Fracture.vertices.col((i - 1) % Fracture.vertices.cols())))
                {
                    estremiTracce.push_back(FirstExtreme);
                    std::cout << FirstExtreme.transpose() << std::endl;
                }
                else if(FractureOperations::isPointOnEdge(SecondExtreme, Fracture.vertices.col(i), Fracture.vertices.col((i - 1) % Fracture.vertices.cols())))
                {
                    estremiTracce.push_back(SecondExtreme);
                    std::cout << SecondExtreme.transpose() << std::endl;
                }
                else
                {
                    bool CurrentCheck = false;
                    Eigen::Vector3d congiungente = Fracture.vertices.col(i) - FirstExtreme;
                    Eigen::Vector3d v= Direction.cross(congiungente);
                    if (v.dot(Fracture.normals)>tol) //controlla > tol oppure >= -tol
                        CurrentCheck = true;
                    if(PreviousCheck != CurrentCheck)
                    {
                        //solve sistem
                        Eigen::MatrixXd A(3,2);
                        A.col(0) = Direction;
                        A.col(1) = Fracture.vertices.col(i) - Fracture.vertices.col((i - 1) % Fracture.vertices.cols());

                        Eigen::Vector3d b = Fracture.vertices.col(i) - FirstExtreme;

                        Eigen::Vector2d paramVert = A.colPivHouseholderQr().solve(b);
                        Eigen::Vector3d Candidate1= Fracture.vertices.col(i) + paramVert[1]*(Fracture.vertices.col((i - 1) % Fracture.vertices.cols())-Fracture.vertices.col(i));

                        Eigen::Vector3d Candidate2 = FirstExtreme + paramVert[0]* Direction;

                        if ( ((Candidate1 - Candidate2)[0] < tol && (Candidate1 - Candidate2)[0] > -tol) &&
                            ((Candidate1 - Candidate2)[1] < tol && (Candidate1 - Candidate2)[1] > -tol) &&
                            ((Candidate1 - Candidate2)[2] < tol && (Candidate1 - Candidate2)[2] > -tol))
                        {
                            estremiTracce.push_back(Candidate1);
                        }
                    }
                }
            }
            //aggiorno FirstExtreme e SecondExtreme con i prolungamenti
            //estremi tracce ha necessariamente due elementi
            std::cout << estremiTracce.size() << std::endl;

            FirstExtreme = estremiTracce[0];
            SecondExtreme = estremiTracce[1];
        }
         //ho prolungato quindi nel caso di tracce passanti ho gli estremi aaltrimenti i prolungamenti
        FirstSide.push_back(FirstExtreme);
        FirstSide.push_back(SecondExtreme);
        SecondSide.push_back(FirstExtreme);
        SecondSide.push_back(SecondExtreme);

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
        }

        //copio std::vector in Eigen::Matrix
        Eigen::MatrixXd FirstSubPolygon(3,FirstSide.size());

        for(unsigned int k = 0; k < FirstSide.size(); k++)
        {
            FirstSubPolygon.col(k) = FirstSide[k];
        }
        std::cout << FirstSubPolygon << std::endl;

        Data::Fract subpolygonuno;
        subpolygonuno.vertices = FirstSubPolygon;
        subpolygonuno.passingTracesId = Fracture.passingTracesId;
        subpolygonuno.notPassingTracesId= Fracture.notPassingTracesId;

        Eigen::MatrixXd SecondSubPolygon(3,SecondSide.size());
        //std::cout << SecondSide.size() << std::endl;
        for(unsigned int k = 0; k < SecondSide.size(); k++)
        {
            SecondSubPolygon.col(k) = SecondSide[k];
        }
        std::cout << SecondSubPolygon << std::endl;

        Data::Fract subpolygondue;
        subpolygondue.vertices = FirstSubPolygon;
        subpolygondue.passingTracesId = Fracture.passingTracesId;
        subpolygondue.notPassingTracesId= Fracture.notPassingTracesId;
        //elimino la traccia considerata
        AllTraces.remove(CurrentTrace.TraceId); //remove elimina tutte le occorrenze... NO PROBLEM: nell nostra lista non ci sono elementi ripetuti
         //richiamo makecut per ogni sotto pologono
        PolygonalMeshLibrary::MakeCuts(subpolygondue, AllTraces, traces, PolygonalMesh, Cell0DId, Cell1DId);
        PolygonalMeshLibrary::MakeCuts(subpolygondue, AllTraces, traces, PolygonalMesh, Cell0DId, Cell1DId);
    }
    return true;
}
}




