#pragma once

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "Utils.hpp"
#include "Fracture.hpp"

using namespace testing;

namespace UnitTesting
{
    TEST(FractureOperationsTEST, fracDistance)
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

    TEST(FractureOperationsTEST, computePlane)
    {

    }

    TEST(FractureOperationsTEST, isTracePassing)
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



    TEST(FractureOperationsTEST, BookCase)
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
