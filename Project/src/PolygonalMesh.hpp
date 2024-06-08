#pragma once

#include <iostream>
#include <vector>
#include <Eigen/Eigen>

namespace PolygonalMeshLibrary{
struct PolygonalMesh{

    unsigned int Num0DsCell;
    std::vector<unsigned int> Id0DsCell;
    std::vector<Eigen::Vector3d> coord0DsCell;

    unsigned int Num1DsCell;
    std::vector<unsigned int> Id1DsCell;
    std::vector<std::array<unsigned int, 2>> coord1DsCell;

    unsigned int Num2DsCell;
    std::vector<std::vector<unsigned int>> Id2DsCell;
    std::vector<std::array<unsigned int, 3>> Cell2DssVertices = {}; ///< Cell2D Vertices indices, size 1 x NumberCell2DVertices[NumberCell2D]
    std::vector<std::array<unsigned int, 3>> Cell2DsEdges = {}; ///< Cell2D Cell1D indices, size 1 x NumberCell2DEdges[NumberCell2D]
};
}

