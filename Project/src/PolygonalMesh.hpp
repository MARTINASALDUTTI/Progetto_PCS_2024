#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

namespace PolygonalMeshLibrary{
struct PolygonalMesh{

    unsigned int Num0DsCell;
    std::vector<unsigned int> Id0DsCell;
    std::vector<Eigen::Vector3d> coord0DsCell; //ci metto dentro un vector per farci operazioni matematiche

    unsigned int Num1DsCell;
    std::vector<unsigned int> Id1DsCell;
    std::vector<std::array<unsigned int, 2>> coord1DsCell; //sono id, non serve fare operazioni matematiche

    unsigned int Num2DsCell;
    std::vector<std::vector<unsigned int>> Id2DsCell; //non so a priori quanti sono e non mi serve accedere per indice (non hanno id)
    std::vector<std::vector<unsigned int>> coord2DsCell;

};
}

