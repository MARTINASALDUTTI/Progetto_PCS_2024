#pragma once

#include <vector>
#include <Eigen/Eigen>


namespace PolygonalMeshLibrary{
struct PolygonalMesh{

    unsigned int Num0DsCell;
    std::unordered_map<unsigned int, Eigen::Vector3d> coord0DsCellMap;
    std::list<Eigen::Vector3d> Vertices_list;

    unsigned int Num1DsCell;
    std::list<std::array<Eigen::Vector3d, 2>> edges_list;
    std::unordered_map<unsigned int, std::array<unsigned int, 2>> Cell1DMap;


    unsigned int Num2DsCell;
    std::list<Eigen::MatrixXd> subpolygons_list;
    std::unordered_map<unsigned int, std::vector<unsigned int>> Cell2DsVertices = {};
    std::unordered_map<unsigned int, std::vector<unsigned int>> Cell2DsEdges = {};
};
}

