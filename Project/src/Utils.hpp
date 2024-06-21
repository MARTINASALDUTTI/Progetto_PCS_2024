#pragma once


#include "Eigen/Eigen"
#include <queue>


#include "Fracture.hpp"
#include "PolygonalMesh.hpp"

namespace Data{
// ImportData(inputFileName, nFracture, Fractures) legge i dati dal file
bool ImportData(const std::string& inputFileName,
                unsigned int& nFracture,
                std::vector<Data::Fract>& Fractures);

bool ExportData(const std::string& outputFileName,
                const std::vector<Data::Trace>& Traces);

bool ExportSecondFile(const std::string& outputFileName,
                      std::vector<Data::Fract>& Fractures,
                      const std::vector<Data::Trace>& Traces);
}

namespace FractureOperations{
//fracDistance return true if the distance between the centres of two fractures is short enough to contain a trace
bool fracDistance(const Eigen::MatrixXd& FirstFracture,
                  const Eigen::MatrixXd& SecondFracture);

//areParallel return true if the fractures lie on two parallel planes
void computePlane(Data::Fract& Fracture);

//findTraces compute traces
bool findTraces(const Data::Fract& FirstFracture,
                const Data::Fract& SecondFracture,
                const Eigen::Vector3d& t,
                Data::Trace& foundTrace);

// findExtreme: compute traces' extreme
void findPosition(const Data::Fract& Fracture,
                  const Eigen::Vector3d& t,
                  const Eigen::Vector3d& P,
                  std::vector<Eigen::VectorXd>& CandidatePoints);

//findExtreme: compute traces' extreme
bool findExtreme(const Eigen::Vector3d& V1,
                 const Eigen::Vector3d& V2,
                 const Eigen::Vector3d& t,
                 const Eigen::Vector3d& P,
                 Eigen::Vector3d& intersection);

bool isPointInPolygon(const Eigen::Vector3d& point,
                      const Eigen::MatrixXd& Polygon,
                      const Eigen::Vector3d& normal);

bool isPointOnEdge(const Eigen::Vector3d& extremePoint,
                   const Eigen::Vector3d& V1,
                   const Eigen::Vector3d& V2);

bool isTracePassing(const Eigen::MatrixXd& fracture,
                    const Eigen::Vector3d& traceStart,
                    const Eigen::Vector3d& traceEnd);

bool bookCase(const Data::Fract& FirstFracture,
              const Data::Fract& SecondFracture,
              Data::Trace& foundTrace);

}

namespace SortLibrary
{
void merge(std::vector<unsigned int>& vecIdTraces,
           const std::vector<Data::Trace>& traces,
           size_t left,
           size_t center,
           size_t right);

void mergesort(std::vector<unsigned int>& vecIdTraces,
               const std::vector<Data::Trace>& traces,
               size_t left,
               size_t right);

void Mergesort(std::vector<unsigned int>& data,
               const std::vector<Data::Trace>& traces);
}

namespace PolygonalMeshLibrary
{
bool MakeCuts(std::list<unsigned int>& AllTraces,
              std::vector<Data::Trace>& traces,
              PolygonalMeshLibrary::PolygonalMesh& PolygonMesh,
              std::queue<Data::Fract>& AllSubPolygons);


bool SolveSystem(const Eigen::Vector3d& Direction,
                 const Eigen::Vector3d& vertice1,
                 const Eigen::Vector3d& vertice2,
                 const Eigen::Vector3d& point,
                 Eigen::Vector3d& Solution);

//void perchè la mesh è da aggiornare
void CreateMesh(PolygonalMeshLibrary::PolygonalMesh& PolygonMesh);

void SavingSubpolygon(const Data::Fract& CurrentPolygon,
                     PolygonalMeshLibrary::PolygonalMesh& PolygonMesh);

bool checking(const Data::Fract& CurrentPolygon,
              Data::Trace& CurrentTrace,
              bool& TraceOnEdge);

bool IsTraceInSubpolygon(const Data::Fract& CurrentPolygon,
                         Data::Trace& CurrentTrace,
                         std::list<unsigned int>& AllTraces,
                         std::vector<Data::Trace>& traces);

bool UpdateTrace(Data::Trace& CurrentTrace,
             std::list<unsigned int>& AllTraces,
             std::vector<Data::Trace>& traces,
             std::vector<Eigen::Vector3d>& estremiTracce);

bool updatetrace(const Eigen::Vector3d V1,
                 const Eigen::Vector3d V2,
                 Data::Trace& CurrentTrace,
                 std::list<unsigned int>& AllTraces,
                 std::vector<Data::Trace>& traces);

bool CalculateArea(const Eigen::MatrixXd& Cell);
}

