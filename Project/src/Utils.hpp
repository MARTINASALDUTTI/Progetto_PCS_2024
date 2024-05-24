#pragma once

#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

#include "Fracture.hpp"

namespace Data{
// ImportData(inputFileName, nFracture, Fractures) legge i dati dal file
bool ImportData(const std::string& inputFileName,
                unsigned int& nFracture,
                std::vector<Data::Fract>& Fractures);
}

namespace FractureOperations{
//fracDistance return true if the distance between the centres of two fractures is short enough to contain a trace
bool fracDistance(Eigen::MatrixXd& FirstFracture,
                  Eigen::MatrixXd& SecondFracture);

//areParallel return true if the fractures lie on two parallel planes
void computePlane(Data::Fract& Fracture);

//findTraces compute traces
void findTraces(const Data::Fract& FirstFracture,
                const Data::Fract& SecondFracture,
                const Eigen::Vector3d& t,
                std::vector<Data::Trace>& Traces);

// findExtreme: compute traces' extreme
void findPosition(const Data::Fract& FirstFracture,
                  const Data::Fract& SecondFracture,
                  const Eigen::Vector3d& t,
                  const Eigen::Vector3d& P,
                  std::vector<Data::Trace>& Traces);

//findExtreme: compute traces' extreme
void findExtreme(const Eigen::Vector3d& V1,
                 const Eigen::Vector3d& V2,
                 const Eigen::Vector3d& t,
                 const Eigen::Vector3d& P,
                 const Data::Fract& SecondFracture,
                 std::vector<Data::Trace>& Traces);

bool isPointInPolygon(const Eigen::Vector3d& point,
                      const Data::Fract& polygon);
}
