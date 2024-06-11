#ifndef __UCD_test_HPP__
#define __UCD_test_HPP__

#include "UCDUtilities.hpp"
#include <gtest/gtest.h>

// ***************************************************************************
TEST(TestUCDUtilities, UCDUtilities_Test0Ds)
{
    const std::string exportFolder = "./";

    Eigen::MatrixXd points(3, 8);
    points.col(0) << 0.8,   0,   0;
    points.col(1) << 1, 0, 0;
    points.col(2) << 1, 1, 0;
    points.col(3) << 0.8, 1, 0;
    points.col(4) << 0.8, 0.5, 0;
    points.col(5) << 0, 1, 0;
    points.col(6) << 0, 0.5, 0;
    points.col(7) << 0, 0, 0;

    Gedim::UCDUtilities exporter;

    exporter.ExportPoints(exportFolder + "Cell0Ds.inp",
                          points);

}
// ***************************************************************************

TEST(TestUCDUtilities, UCDUtilities_Test1Ds)
{
    const std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;

    Eigen::MatrixXd points(3, 8);
    points.col(0) << 0.8,   0,   0;
    points.col(1) << 1, 0, 0;
    points.col(2) << 1, 1, 0;
    points.col(3) << 0.8, 1, 0;
    points.col(4) << 0.8, 0.5, 0;
    points.col(5) << 0, 1, 0;
    points.col(6) << 0, 0.5, 0;
    points.col(7) << 0, 0, 0;

    Eigen::MatrixXi edges(2, 11);
    edges.col(0) << 0, 3;
    edges.col(1) << 1, 0;
    edges.col(2) << 2, 1;
    edges.col(3) << 3, 2;
    edges.col(4) << 4, 6;
    edges.col(5) << 5, 4;
    edges.col(6) << 3, 5;
    edges.col(7) << 6, 3;
    edges.col(8) << 0, 4;
    edges.col(9) << 7, 0;
    edges.col(10) << 6, 7;

    exporter.ExportSegments(exportFolder + "Cell1Ds.inp",
                            points,
                            edges);

}
// ***************************************************************************
TEST(TestUCDUtilities, UCDUtilities_Test2Ds)
{
    const std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;

    Eigen::MatrixXd points(3, 8);
    points.col(0) << 0.8,   0,   0;
    points.col(1) << 1, 0, 0;
    points.col(2) << 1, 1, 0;
    points.col(3) << 0.8, 1, 0;
    points.col(4) << 0.8, 0.5, 0;
    points.col(5) << 0, 1, 0;
    points.col(6) << 0, 0.5, 0;
    points.col(7) << 0, 0, 0;

    std::vector<std::vector<unsigned int>> polygons_vertices = {{0, 1, 2, 3},
                                                                {4, 5, 3, 6},
                                                                {4, 0, 7, 6}};

    exporter.ExportPolygons(exportFolder + "Cell2Ds.inp",
                            points,
                            polygons_vertices);

}
// ***************************************************************************

#endif // __UCD_test_HPP__
