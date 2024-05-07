#pragma once

#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

// ImportData(inputFileName, nFracture, Fractures) legge i dati dal file
bool ImportData(const std::string& inputFileName,
                unsigned int& nFracture,
                std::map<unsigned int, Eigen::MatrixXd> Fractures);