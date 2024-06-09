#pragma once

#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

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
bool fracDistance(Eigen::MatrixXd& FirstFracture,
                  Eigen::MatrixXd& SecondFracture);

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
                      const Eigen::MatrixXd& Fracture);

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
bool MakeCuts(const Data::Fract& Fracture,
              std::list<unsigned int>& AllTraces,
              const std::vector<Data::Trace>& traces,
              PolygonalMeshLibrary::PolygonalMesh& PolygonalMesh,
              unsigned int& Cell0DId,
              unsigned int& Cell1DId);

void MakeMesh(const Data::Fract& Fracture,
              PolygonalMeshLibrary::PolygonalMesh& PolygonalMesh,
              //sicuramente non modo migliroe
              unsigned int& Cell0DId,
              unsigned int& Cell1DId);
}
