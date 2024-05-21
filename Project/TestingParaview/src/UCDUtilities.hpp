
#ifndef __TEST_EXPORT_POLYGONS_TO_PARAVIEW_H
#define __TEST_EXPORT_POLYGONS_TO_PARAVIEW_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <iostream>
#include <vector>

#include "export_polygons_to_paraview.hpp"
#include "IOUtilities.hpp"

using namespace testing;

namespace UnitTesting
{
// ***************************************************************************
TEST(TestExportPolygons, TestExportPolygonsToParaview)
{
    const std::string exportFolder = "./Export";
    Gedim::Output::CreateFolder(exportFolder);

    Eigen::MatrixXd points(3, 4);
    points.col(0)<< 0.0000000000000000e+00; 1.0000000000000000e+00; 1.0000000000000000e+00; 0.0000000000000000e+00
;
    points.col(1)<< 1.0, 0.0, 0.0;
    points.col(2)<< 0.0, 1.0, 0.0;
    points.col(3)<< 0.0, 0.0, 1.0;

    std::vector<std::vector<unsigned int>> polygons_vertices = {{0, 1, 2, 3},
                                                                {2, 1, 0}};

    const std::string file_path = exportFolder + "/polygons.inp";

    Functions::export_polygons_to_paraview(points,
                                           polygons_vertices,
                                           file_path);

}
// ***************************************************************************
}

#endif // __TEST_EXPORT_POLYGONS_TO_PARAVIEW_H
