#pragma once

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "Utils.hpp"
#include "Fracture.hpp"

constexpr double tol = 10e-10;

using namespace testing;

namespace FractureOperations
{
    TEST(FractureOperationsTEST, testfracDistance)
    {
        Eigen::MatrixXd A(3,4);
        A << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

        Eigen::MatrixXd B(3,4);
        B << 8.0000000000000004e-01, 8.0000000000000004e-01, 8.0000000000000004e-01, 8.0000000000000004e-01,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            -1.0000000000000001e-01, 2.9999999999999999e-01, 2.9999999999999999e-01, -1.0000000000000001e-01;

        EXPECT_TRUE(FractureOperations::fracDistance(A,B));
    }

    TEST(FractureOperationsTEST, testcomputePlane)
    {
        Data::Fract TestFract;
        TestFract.FractId = 0;
        Eigen::MatrixXd A(3,4);
        A << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
        TestFract.vertices = A;
        FractureOperations::computePlane(TestFract);
        Eigen::Vector3d normal;
        normal << 0, 0, 1;
        double knownTerm = 0;
        ASSERT_TRUE(TestFract.normals == normal);
        ASSERT_TRUE(TestFract.d == knownTerm);
    }

    TEST(FractureOperationsTEST, testfindTraces)
    {
        Data::Fract TestFract1;
        Data::Fract TestFract2;
        Eigen::MatrixXd A(3,4);
        A << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
        Eigen::MatrixXd B(3,4);
        B << 8.0000000000000004e-01, 8.0000000000000004e-01, 8.0000000000000004e-01, 8.0000000000000004e-01,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            -1.0000000000000001e-01, 2.9999999999999999e-01, 2.9999999999999999e-01, -1.0000000000000001e-01;
        TestFract1.FractId = 0;
        TestFract2.FractId = 1;
        TestFract1.vertices = A;
        TestFract2.vertices = B;
        Eigen::Vector3d normal1;
        normal1 << 0, 0, 1;
        TestFract1.normals = normal1;
        Eigen::Vector3d normal2;
        normal2 << -1, 0, 0;
        TestFract2.normals = normal2;
        double knownTerm1 = 0;
        TestFract1.d = knownTerm1;
        double knownTerm2 = -0.8;
        TestFract2.d = knownTerm2;
        Eigen::Vector3d t;
        t << 0, -1, 0;
        Data::Trace testTrace;
        ASSERT_TRUE(FractureOperations::findTraces(TestFract1, TestFract2, t, testTrace));
        std::vector<Eigen::Vector3d> ExtremesTest;
        Eigen::Vector3d a, b;
        a << 0.8, 0, 0;
        b << 0.8, 1, 0;
        ExtremesTest.push_back(a);
        ExtremesTest.push_back(b);
        ASSERT_TRUE(testTrace.ExtremesCoord[0] == ExtremesTest[0] && testTrace.ExtremesCoord[1] == ExtremesTest[1]);
        ASSERT_TRUE(testTrace.Tips[0] == false);
        ASSERT_TRUE(testTrace.Tips[1] == false);
        ASSERT_TRUE(testTrace.length - 1 <= tol && testTrace.length - 1  >= -tol);
    }

    TEST(FractureOperationsTEST, testisTracePassing)
    {
        Eigen::MatrixXd fracture(3,4);
        fracture << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

        Eigen::Vector3d traceStart;
        traceStart << 0.8, 0, 0;

        Eigen::Vector3d traceEnd;
        traceEnd << 0.8, 1, 0;

        ASSERT_TRUE(FractureOperations::isTracePassing(fracture, traceStart, traceEnd));

        Eigen::Vector3d SecondtraceStart;
        SecondtraceStart << 0.5, 0, 0;

        Eigen::Vector3d SecondtraceEnd;
        SecondtraceEnd << 0.316184, 0.5, 0;

        ASSERT_FALSE(FractureOperations::isTracePassing(fracture, SecondtraceStart, SecondtraceEnd));
    }

    TEST(FractureOperationsTEST, testisPointInPolygon)
    {
        Eigen::Vector3d Firstpoint;
        Firstpoint << 0, 0.5, 0;
        Eigen::MatrixXd Fracture(3,4);
        Fracture << 8.0000000000000004e-01, 8.0000000000000004e-01, 8.0000000000000004e-01, 8.0000000000000004e-01,
                    0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
                   -1.0000000000000001e-01, 2.9999999999999999e-01, 2.9999999999999999e-01, -1.0000000000000001e-01;
        Eigen::Vector3d normal;
        normal << -1, 0, 0;
        ASSERT_TRUE(FractureOperations::isPointInPolygon(Firstpoint, Fracture, normal));

        Eigen::Vector3d Secondpoint;
        Secondpoint << 1, 1, 1;
        ASSERT_FALSE(FractureOperations::isPointInPolygon(Secondpoint, Fracture, normal));

    }

    TEST(FractureOperationsTEST, testisPointOnEdge)
    {
        Eigen::Vector3d FirstExtremePoint;
        FirstExtremePoint << 0.8, 0, 0;
        Eigen::Vector3d V1;
        V1 << 0, 0, 0;
        Eigen::Vector3d V2;
        V2 << 1, 0, 0;
        ASSERT_TRUE(FractureOperations::isPointOnEdge(FirstExtremePoint, V1, V2));

        Eigen::Vector3d SecondExtremePoint;
        SecondExtremePoint << -0.8, 0, 0;
        ASSERT_FALSE(FractureOperations::isPointOnEdge(SecondExtremePoint, V1, V2));

    }

    TEST(FractureOperationsTEST, testBookCase)
    {
        Data::Fract FirstFracture;
        Eigen::MatrixXd A(3,4);
        A << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
        FirstFracture.vertices = A;
        FirstFracture.FractId = 0;

        Data::Fract SecondFracture;
        Eigen::MatrixXd B(3,4);
        B << 0.0000000000000000e+00, 1.0000000000000000e+00, -1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, -1.0000000000000000e+00, -1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

        SecondFracture.vertices = B;
        SecondFracture.FractId = 1;

        Data::Trace foundTrace;
        ASSERT_FALSE(FractureOperations::bookCase(FirstFracture, SecondFracture,foundTrace));

    }
}

namespace SortLibrary
{
    TEST(SortLibraryTEST, testMergeSort)
    {
        std::vector<Data::Trace> vector_testMergeSort(6);
        Data::Trace trace0;
        trace0.TraceId = 0;
        trace0.length = 3.0;
        Data::Trace trace1;
        trace1.TraceId = 1;
        trace1.length = 1.0;
        Data::Trace trace2;
        trace2.TraceId = 2;
        trace2.length = 0.5;
        Data::Trace trace3;
        trace3.TraceId = 3;
        trace3.length = 8.0;
        Data::Trace trace4;
        trace4.TraceId = 4;
        trace4.length = 5.0;
        Data::Trace trace5;
        trace5.TraceId = 5;
        trace5.length = 2.0;
        vector_testMergeSort[0] = trace0;
        vector_testMergeSort[1] = trace1;
        vector_testMergeSort[2] = trace2;
        vector_testMergeSort[3] = trace3;
        vector_testMergeSort[4] = trace4;
        vector_testMergeSort[5] = trace5;

        std::vector<unsigned int> vectorID_testMergeSort = {0, 1, 2, 3, 4, 5};

        SortLibrary::Mergesort(vectorID_testMergeSort, vector_testMergeSort);

        EXPECT_EQ(vectorID_testMergeSort[0], 3);
        EXPECT_EQ(vectorID_testMergeSort[1], 4);
        EXPECT_EQ(vectorID_testMergeSort[2], 0);
        EXPECT_EQ(vectorID_testMergeSort[3], 5);
        EXPECT_EQ(vectorID_testMergeSort[4], 1);
        EXPECT_EQ(vectorID_testMergeSort[5], 2);
    }
}
