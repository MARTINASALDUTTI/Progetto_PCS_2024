#pragma once

#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

#include "Fracture.hpp"

// ImportData(inputFileName, nFracture, Fractures) legge i dati dal file
bool ImportData(const std::string& inputFileName,
                unsigned int& nFracture,
                std::vector<Eigen::MatrixXd>& Fractures);

//fracDistance return true if the distance between the centres of two fractures is short enough to contain a trace
bool fracDistance(Eigen::MatrixXd& FirstFracture,
                  Eigen::MatrixXd& SecondFracture);

//areParallel return true if the fractures lie on two parallel planes
bool areParallel(Eigen::MatrixXd& FirstFracture,
                 Eigen::MatrixXd& SecondFracture);

//findTraces compute traces
void findTraces(Eigen::MatrixXd& FirstFracture,
                Eigen::MatrixXd& SecondFracture,
                std::vector<Trace>& Traces);
