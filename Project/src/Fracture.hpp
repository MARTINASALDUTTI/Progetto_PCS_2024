#pragma once

#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

namespace Data{
struct Fract {
    unsigned int FractId;
    Eigen::MatrixXd vertices;
    Eigen::Vector3d normals;
    double d; //known term of the fracture
    std::vector<unsigned int> passingTracesId;  // Id of passing traces in the fracture
    std::vector<unsigned int> notPassingTracesId;    // Id of non passing traces in the fracture
};

struct Trace{
    unsigned int TraceId;
    std::array<unsigned int, 2> FractureIds={};
    std::array<Eigen::Vector3d, 2> ExtremesCoord={}; //array of coordinates of traces' extremes
    std::array<bool,2> Tips={}; //false for passing traces and true for not passing traces
    double length;
    //std::array<bool,2> passThrough;
};
}
