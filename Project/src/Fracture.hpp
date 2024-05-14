#pragma once

#include <iostream>
#include <sstream>
#include "Eigen/Eigen"

struct Trace{
    unsigned int TraceId;
    std::array<unsigned int, 2> FractureIds={};
    std::array<Eigen::Vector3d, 2> ExtremesCoord={}; //array of coordinates of traces' extremes
    std::array<bool,2> Tips={}; //false for passing traces and true for not passing traces
    double length;
    //std::array<bool,2> passThrough;
};
