#pragma once

#include "Eigen/Eigen"
#include "UCDUtilities.hpp"

namespace Functions
{
// ***************************************************************************
void export_polygons_to_paraview(const Eigen::MatrixXd& points,
                                 const std::vector<std::vector<unsigned int>>& polygons_vertices,
                                 const std::string& file_path)
{
    Gedim::UCDUtilities exporter;

    exporter.ExportPolygons(file_path,
                            points,
                            polygons_vertices);

}
// ***************************************************************************
}
