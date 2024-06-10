#ifndef __UCD_test_HPP__
#define __UCD_test_HPP__

#include "UCDUtilities.hpp"
#include <gtest/gtest.h>

// ***************************************************************************
TEST(TestUCDUtilities, UCDUtilities_Test0Ds)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;

    exporter.ExportPoints(exportFolder + "/Cell0Ds.inp",
                          points);

}
// ***************************************************************************
TEST(TestUCDUtilities, UCDUtilities_Test1Ds)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;


    exporter.ExportSegments(exportFolder + "/Cell1Ds.inp",
                            points,
                            edges);
}
// ***************************************************************************
TEST(TestUCDUtilities, UCDUtilities_Test2Ds)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;

    exporter.ExportPolygons(exportFolder + "/Cell2Ds.inp",
                            points,
                            polygons);
}
// ***************************************************************************

#endif // __UCD_test_HPP__
