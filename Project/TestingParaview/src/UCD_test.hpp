#ifndef __UCD_test_HPP__
#define __UCD_test_HPP__

#include "UCDUtilities.hpp"
#include <gtest/gtest.h>

// ***************************************************************************
TEST(TestUCDUtilities, UCDUtilities_Test0Ds)
{
    const std::string exportFolder = "./";

    /*
    Eigen::MatrixXd points(3, 8);
    points.col(0) << 0.8,   0,   0;
    points.col(1) << 1, 0, 0;
    points.col(2) << 1, 1, 0;
    points.col(3) << 0.8, 1, 0;
    points.col(4) << 0, 0.5, 0;
    points.col(5) << 0.8, 0.5, 0;
    points.col(6) << 0, 1, 0;
    points.col(7) << 0, 0, 0;
    */

    Eigen::MatrixXd points(3, 10);
    points.row(0) << 6.7949650570084286e-01, 2.0959413133064569e-01, 7.7027229455623514e-02, 5.4692960382582068e-01, 0.554138, 0.422841, 0.160363,  0.424368, 0.407083,  0.589994;
    points.row(1) << 5.1566886122886846e-01, 9.9389350435296486e-01, 8.6363358811981283e-01, 3.8540894499571632e-01, 0.392492, 0.77687,  0.945519, 0.544407, 0.792907,  0.441441;
    points.row(2) << 1.9054542365205804e-01, 1.9054542365205804e-01, 5.9704177318758545e-01, 5.9704177318758545e-01, 0.574938, 0.190545,  0.341506, 0.544514, 0.190545,  0.443965;

    Gedim::UCDUtilities exporter;

    exporter.ExportPoints(exportFolder + "Cell0Ds.inp",
                          points);

}
// ***************************************************************************

TEST(TestUCDUtilities, UCDUtilities_Test1Ds)
{
    const std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;

    /*
    Eigen::MatrixXd points(3, 8);
    points.col(0) << 0.8,   0,   0;
    points.col(1) << 1, 0, 0;
    points.col(2) << 1, 1, 0;
    points.col(3) << 0.8, 1, 0;
    points.col(4) << 0, 0.5, 0;
    points.col(5) << 0.8, 0.5, 0;
    points.col(6) << 0, 1, 0;
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
    edges.col(8) << 7, 4;
    edges.col(9) << 0, 7;
    edges.col(10) << 5, 0;
*/
    Eigen::MatrixXd points(3, 10);
    points.row(0) << 6.7949650570084286e-01, 2.0959413133064569e-01, 7.7027229455623514e-02, 5.4692960382582068e-01, 0.554138, 0.422841, 0.160363,  0.424368, 0.407083,  0.589994;
    points.row(1) << 5.1566886122886846e-01, 9.9389350435296486e-01, 8.6363358811981283e-01, 3.8540894499571632e-01, 0.392492, 0.77687,  0.945519, 0.544407, 0.792907,  0.441441;
    points.row(2) << 1.9054542365205804e-01, 1.9054542365205804e-01, 5.9704177318758545e-01, 5.9704177318758545e-01, 0.574938, 0.190545,  0.341506, 0.544514, 0.190545,  0.443965;

    Eigen::MatrixXi edges(2,7);
    edges.col(0) << 0, 1;
    edges.col(1) << 1, 2;
    edges.col(2) << 2, 3;
    edges.col(3) << 3, 0;
    edges.col(4) << 4, 5;
    edges.col(5) << 6, 7;
    edges.col(6) << 8, 9;

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
    points.col(4) << 0, 0.5, 0;
    points.col(5) << 0.8, 0.5, 0;
    points.col(6) << 0, 1, 0;
    points.col(7) << 0, 0, 0;

    std::vector<std::vector<unsigned int>> polygons_vertices = {{0, 1, 2, 3},
                                                                {4, 5, 3, 6},
                                                                {7, 0, 5, 4}};

    exporter.ExportPolygons(exportFolder + "Cell2Ds.inp",
                            points,
                            polygons_vertices);

}
// ***************************************************************************

#endif // __UCD_test_HPP__
