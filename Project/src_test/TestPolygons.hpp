#pragma once

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>

#include "Utils.hpp"
#include "Fracture.hpp"

using namespace testing;

namespace UnitTesting
{

    TEST(FractureOperationTEST, fracDistance)
    {
        Eigen::MatrixXd A;
        A << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

        Eigen::MatrixXd B;
        B << 8.0000000000000004e-01, 8.0000000000000004e-01, 8.0000000000000004e-01, 8.0000000000000004e-01,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            -1.0000000000000001e-01, 2.9999999999999999e-01, 2.9999999999999999e-01, -1.0000000000000001e-01;

        EXPECT_TRUE(FractureOperations::fracDistance(A,B));
    }

    TEST(FractureOperationTEST, computePlane)
    {

    }


    TEST(TestFractureOperations, BookCase)
    {
        Data::Fract FirstFracture;
        Eigen::MatrixXd A;
        A << 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 1.0000000000000000e+00, 1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;
        FirstFracture.vertices = A;


        Data::Fract SecondFracture;
        Eigen::MatrixXd B;
        B << 0.0000000000000000e+00, 1.0000000000000000e+00, -1.0000000000000000e+00, 0.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, -1.0000000000000000e+00, -1.0000000000000000e+00,
            0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00;

        SecondFracture.vertices = B;

        Data::Trace foundTrace;

        FractureOperations::bookCase(FirstFracture, SecondFracture,foundTrace);
    }
}
