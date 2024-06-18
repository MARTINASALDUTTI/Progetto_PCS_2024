#ifndef __UCD_test_HPP__
#define __UCD_test_HPP__

#include "UCDUtilities.hpp"
#include <gtest/gtest.h>

// ***************************************************************************
TEST(TestUCDUtilities, Test0Ds)
{
    const std::string exportFolder = "./";


    Eigen::MatrixXd points(3, 26);
    points.col(0) << 0.594717, 0.432365, 0.450509;
    points.col(1) << 0.450525, 0.709432, 0.250734;
    points.col(2) << 0.446922, 0.706373, 0.261045;
    points.col(3) << 0.554138, 0.392492, 0.574938;
    points.col(4) << 0.436999, 0.735424, 0.231993;
    points.col(5) << 0.679497, 0.515669, 0.190545;
    points.col(6) << 0.479352, 0.719358, 0.190545;
    points.col(7) << 0.433893, 0.744516, 0.2229;
    points.col(8) << 0.422841, 0.77687, 0.190545;
    points.col(9) << 0.364069, 0.636022, 0.498147;
    points.col(10) << 0.492678, 0.440621, 0.597042;
    points.col(11) << 0.54693, 0.385409, 0.597042;
    points.col(12) << 0.431609, 0.74578, 0.224526;
    points.col(13) << 0.247976, 0.689657, 0.597042;
    points.col(14) << 0.241961, 0.821544, 0.404251;
    points.col(15) << 0.160363, 0.945519, 0.341506;
    points.col(16) << 0.0770272, 0.863634, 0.597042;
    points.col(17) << 0.329511, 0.606678, 0.597042;
    points.col(18) << 0.4272, 0.507259, 0.597042;
    points.col(19) << 0.352669, 0.626342, 0.530771;
    points.col(20) << 0.241961, 0.821544, 0.404251;
    points.col(21) << 0.235293, 0.967739, 0.190545;
    points.col(22) << 0.209594, 0.993894, 0.190545;
    points.col(23) << 0.407083, 0.792907, 0.190545;
    points.col(24) << 0.433893, 0.744516, 0.2229;
    points.col(25) << 0.240595, 0.851491, 0.360475;


    Gedim::UCDUtilities exporter;

    exporter.ExportPoints(exportFolder + "Cell0Ds.inp",
                          points);

}
// ***************************************************************************

TEST(TestUCDUtilities, Test1Ds)
{
    const std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;

    Eigen::MatrixXd points(3, 4);
    points << 6.7949650570084286e-01,  2.0959413133064569e-01,  7.7027229455623514e-02,  5.4692960382582068e-01,
              5.1566886122886846e-01,  9.9389350435296486e-01,  8.6363358811981283e-01,  3.8540894499571632e-01,
              1.9054542365205804e-01,  1.9054542365205804e-01,  5.9704177318758545e-01,  5.9704177318758545e-01;


    Eigen::MatrixXi edges(2, 4);
    edges.col(0) << 0, 1;
    edges.col(1) << 1, 2;
    edges.col(2) << 2, 3;
    edges.col(3) << 3, 0;

    exporter.ExportSegments(exportFolder + "Fracture.inp",
                        points,
                        edges);

}
// ***************************************************************************

TEST(TestUCDUtilities, Test1DsTrace)
{
    const std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    Eigen::MatrixXd points(3, 14);

    points.col(0) << 0.554138, 0.392492, 0.574938;
    points.col(1) <<    0.422841, 0.77687, 0.190545;
    points.col(2) <<    0.160363, 0.945519, 0.341506;
    points.col(3) <<    0.424368, 0.544407, 0.544514;
    points.col(4) <<    0.407083, 0.792907, 0.190545;
    points.col(5) <<    0.589994, 0.441441, 0.443965;
    points.col(6) <<    0.329511, 0.606678, 0.597042;
    points.col(7) <<    0.449954, 0.708947, 0.252369;
    points.col(8) <<    0.247976, 0.689657, 0.597042;
    points.col(9) <<    0.238807, 0.890697, 0.303164;
    points.col(10) <<    0.479352, 0.719358, 0.190545;
    points.col(11) <<    0.322542, 0.80614, 0.302152;
    points.col(12) <<    0.4272, 0.507259, 0.597042;
    points.col(13) <<    0.425099, 0.510615, 0.595174;

    Eigen::MatrixXi edges(2, 7);
    edges.col(0) << 0, 1;
    edges.col(1) << 2, 3;
    edges.col(2) << 4, 5;
    edges.col(3) << 6, 7;
    edges.col(4) << 8, 9;
    edges.col(5) << 10, 11;
    edges.col(6) << 12, 13;

    exporter.ExportSegments(exportFolder + "Trace.inp",
                            points,
                            edges);

}

TEST(TestUCDUtilities, Test1DsCell)
{
    const std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;

    Eigen::MatrixXd points(3, 26);
    points.col(0) << 0.594717, 0.432365, 0.450509;
    points.col(1) << 0.450525, 0.709432, 0.250734;
    points.col(2) << 0.446922, 0.706373, 0.261045;
    points.col(3) << 0.554138, 0.392492, 0.574938;
    points.col(4) << 0.436999, 0.735424, 0.231993;
    points.col(5) << 0.679497, 0.515669, 0.190545;
    points.col(6) << 0.479352, 0.719358, 0.190545;
    points.col(7) << 0.433893, 0.744516, 0.2229;
    points.col(8) << 0.422841, 0.77687, 0.190545;
    points.col(9) << 0.364069, 0.636022, 0.498147;
    points.col(10) << 0.492678, 0.440621, 0.597042;
    points.col(11) << 0.54693, 0.385409, 0.597042;
    points.col(12) << 0.431609, 0.74578, 0.224526;
    points.col(13) << 0.247976, 0.689657, 0.597042;
    points.col(14) << 0.241961, 0.821544, 0.404251;
    points.col(15) << 0.160363, 0.945519, 0.341506;
    points.col(16) << 0.0770272, 0.863634, 0.597042;
    points.col(17) << 0.329511, 0.606678, 0.597042;
    points.col(18) << 0.4272, 0.507259, 0.597042;
    points.col(19) << 0.352669, 0.626342, 0.530771;
    points.col(20) << 0.241961, 0.821544, 0.404251;
    points.col(21) << 0.235293, 0.967739, 0.190545;
    points.col(22) << 0.209594, 0.993894, 0.190545;
    points.col(23) << 0.407083, 0.792907, 0.190545;
    points.col(24) << 0.433893, 0.744516, 0.2229;
    points.col(25) << 0.240595, 0.851491, 0.360475;

    Eigen::MatrixXi edges(2, 44);

    edges.col(0) << 0, 1;
    edges.col(1) << 1, 2;
    edges.col(2) << 2, 3;
    edges.col(3) << 3, 0;
    edges.col(4) << 1, 4;
    edges.col(5) << 4, 2;
    edges.col(6) << 5, 6;
    edges.col(7) << 6, 7;
    edges.col(8) << 7, 4;
    edges.col(9) << 4,0;
    edges.col(10) << 0, 5;
    edges.col(11) << 6, 8;
    edges.col(12) << 8 , 7;
    edges.col(13) << 2, 9;
    edges.col(14) << 9, 10;
    edges.col(15) << 10, 11;
    edges.col(16) << 11, 3;
    edges.col(17) << 7, 12;
    edges.col(18) << 12, 4;
    edges.col(19) << 13, 14;
    edges.col(20) << 14, 15;
    edges.col(21) << 15, 16;
    edges.col(22) << 16, 13;
    edges.col(23) << 17, 9;
    edges.col(24) << 9, 14;
    edges.col(25) << 13, 17;
    edges.col(26) << 18, 19;
    edges.col(27) << 19, 17;
    edges.col(28) << 17, 18;
    edges.col(29) << 9, 19;
    edges.col(30) << 18, 10;
    edges.col(31) << 20, 21;
    edges.col(32) << 21, 22;
    edges.col(33) << 22, 15;
    edges.col(34) << 15, 20;
    edges.col(35) << 12, 23;
    edges.col(36) << 24, 12;
    edges.col(37) << 8, 24;
    edges.col(38) << 23, 8;
    edges.col(39) << 12, 25;
    edges.col(40) << 20, 25;
    edges.col(41) << 9, 20;
    edges.col(42) << 23, 21;
    edges.col(43) << 25, 21;

    exporter.ExportSegments(exportFolder + "Cell1Ds.inp",
                            points,
                            edges);
}

// ***************************************************************************

#endif // __UCD_test_HPP__
